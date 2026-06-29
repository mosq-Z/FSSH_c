#include "fssh_thread_pool.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <complex.h>
#include <windows.h>

static int get_cpu_core_count(void)
{
    SYSTEM_INFO sys_info;
    GetSystemInfo(&sys_info);
    return (int)sys_info.dwNumberOfProcessors;
}

// 工作线程主函数
static void* worker_thread_func(void* arg)
{
    ThreadPool* pool = (ThreadPool*)arg;
    // 每个线程设置独立随机种子，解决多线程rand()线程不安全问题
    unsigned int seed = (unsigned int)(time(NULL) + pthread_self());
    srand(seed);

    while (true)
    {
        pthread_mutex_lock(&pool->mutex);

        // 无任务则等待条件变量
        while (pool->task_head == NULL && !pool->pool_shutdown)
        {
            pthread_cond_wait(&pool->cond, &pool->mutex);
        }

        // 线程池关闭，退出工作线程
        if (pool->pool_shutdown)
        {
            pthread_mutex_unlock(&pool->mutex);
            pthread_exit(NULL);
        }

        // 取出队首任务
        TaskNode* node = pool->task_head;
        pool->task_head = pool->task_head->next;
        pthread_mutex_unlock(&pool->mutex);

        // 执行FSSH单条轨迹
        TrajTask task = node->task;
        free(node);

        TrajContext traj;
        traj_init(&traj, task.k);
        traj_run_full(&traj);

        int traj_res[2];
        traj_get_result(&traj, traj_res);

        // 写入结果数组（加互斥锁保证线程安全）
        pthread_mutex_lock(&pool->res_mutex);
        pool->result_arr[pool->finished_cnt].traj_idx = task.traj_idx;
        pool->result_arr[pool->finished_cnt].surface_flag = traj_res[0];
        pool->result_arr[pool->finished_cnt].trans_flag = traj_res[1];
        pool->finished_cnt++;
        pthread_mutex_unlock(&pool->res_mutex);
    }
    return NULL;
}

// 创建线程池
ThreadPool* thread_pool_create(void)
{
    ThreadPool* pool = (ThreadPool*)malloc(sizeof(ThreadPool));
    if (!pool) return NULL;

    // 初始化锁与条件变量
    pthread_mutex_init(&pool->mutex, NULL);
    pthread_cond_init(&pool->cond, NULL);
    pthread_mutex_init(&pool->res_mutex, NULL);

    pool->task_head = NULL;
    pool->pool_shutdown = false;
    pool->finished_cnt = 0;

    // 自动获取CPU物理核心数
    pool->thread_num = get_cpu_core_count();
    pool->workers = (pthread_t*)malloc(sizeof(pthread_t) * pool->thread_num);

    // 预分配结果数组空间
    long total_traj = g_fssh_cfg.N_Single_K_Traj;
    pool->result_arr = (TrajTaskResult*)malloc(sizeof(TrajTaskResult) * total_traj);

    // 创建所有工作线程
    for (int i = 0; i < pool->thread_num; i++)
    {
        pthread_create(&pool->workers[i], NULL, worker_thread_func, pool);
    }

    printf("线程池创建成功，工作线程数：%d（CPU核心数）\n", pool->thread_num);
    return pool;
}

// 提交任务到任务队列
int thread_pool_submit_task(ThreadPool* pool, TrajTask task)
{
    if (!pool || pool->pool_shutdown) return -1;

    TaskNode* new_node = (TaskNode*)malloc(sizeof(TaskNode));
    new_node->task = task;
    new_node->next = NULL;

    pthread_mutex_lock(&pool->mutex);
    if (pool->task_head == NULL)
    {
        pool->task_head = new_node;
    }
    else
    {
        TaskNode* cur = pool->task_head;
        while (cur->next != NULL)
            cur = cur->next;
        cur->next = new_node;
    }
    pthread_cond_signal(&pool->cond);
    pthread_mutex_unlock(&pool->mutex);

    return 0;
}

// 阻塞等待所有任务完成
void thread_pool_wait_all(ThreadPool* pool)
{
    while (true)
    {
        pthread_mutex_lock(&pool->mutex);
        bool empty = (pool->task_head == NULL);
        pthread_mutex_unlock(&pool->mutex);

        if (empty && pool->finished_cnt == g_fssh_cfg.N_Single_K_Traj)
            break;
        usleep(10000);
    }
}

// 销毁线程池
void thread_pool_destroy(ThreadPool* pool)
{
    if (!pool) return;

    pthread_mutex_lock(&pool->mutex);
    pool->pool_shutdown = true;
    pthread_cond_broadcast(&pool->cond);
    pthread_mutex_unlock(&pool->mutex);

    // 等待所有线程退出
    for (int i = 0; i < pool->thread_num; i++)
    {
        pthread_join(pool->workers[i], NULL);
    }

    // 释放剩余任务节点
    TaskNode* cur = pool->task_head;
    while (cur)
    {
        TaskNode* tmp = cur;
        cur = cur->next;
        free(tmp);
    }

    // 释放资源
    free(pool->workers);
    free(pool->result_arr);
    pthread_mutex_destroy(&pool->mutex);
    pthread_cond_destroy(&pool->cond);
    pthread_mutex_destroy(&pool->res_mutex);
    free(pool);
}

// 单动量多核并行批量模拟
SingleKResult batch_run_single_k_parallel(double k)
{
    SingleKResult res;
    res.k = k;
    long total = g_fssh_cfg.N_Single_K_Traj;

    // 创建线程池
    ThreadPool* pool = thread_pool_create();

    // 批量提交所有轨迹任务
    for (long i = 0; i < total; i++)
    {
        TrajTask task = {
            .k = k,
            .traj_idx = (int)i
        };
        thread_pool_submit_task(pool, task);
    }

    // 等待全部轨迹计算完成
    thread_pool_wait_all(pool);

    // 统计四类结果
    long count_TL = 0, count_RL = 0;
    long count_TH = 0, count_RH = 0;
    long count_other = 0;

    for (long i = 0; i < total; i++)
    {
        TrajTaskResult tr = pool->result_arr[i];
        int surface = tr.surface_flag;
        int trans = tr.trans_flag;

        if (surface == 0)
        {
            if (trans == 1) count_TL++;
            else count_RL++;
        }
        else
        {
            if (trans == 1) count_TH++;
            else count_RH++;
        }
    }

    // 归一化概率
    double total_t = (double)total;
    res.TL = (double)count_TL / total_t;
    res.RL = (double)count_RL / total_t;
    res.TH = (double)count_TH / total_t;
    res.RH = (double)count_RH / total_t;
    res.other = (double)count_other / total_t;

    // 销毁线程池释放资源
    thread_pool_destroy(pool);
    return res;
}