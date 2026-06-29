#ifndef FSSH_THREAD_POOL_H
#define FSSH_THREAD_POOL_H

#include "fssh_batch.h"
#include <pthread.h>
#include <stdbool.h>

// 单条轨迹任务：仅保存当前轨迹需要的参数
typedef struct {
    double k;
    int traj_idx;
} TrajTask;

// 单条轨迹运行完返回的结果
typedef struct {
    int traj_idx;
    int surface_flag;   // res[0] 0低能/1高能
    int trans_flag;     // res[1] 1透射/0反射
} TrajTaskResult;

// 任务队列节点
typedef struct TaskNode {
    TrajTask task;
    struct TaskNode* next;
} TaskNode;

// 线程池结构体
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;

    TaskNode* task_head;
    bool pool_shutdown;

    int thread_num;
    pthread_t* workers;

    // 结果存储 + 结果互斥锁
    TrajTaskResult* result_arr;
    pthread_mutex_t res_mutex;
    long finished_cnt;
} ThreadPool;

/**
 * @brief 创建线程池，自动获取CPU核心数创建工作线程
 * @return 线程池指针
 */
ThreadPool* thread_pool_create(void);

/**
 * @brief 向线程池提交一条轨迹任务
 */
int thread_pool_submit_task(ThreadPool* pool, TrajTask task);

/**
 * @brief 等待所有任务执行完成，阻塞主线程
 */
void thread_pool_wait_all(ThreadPool* pool);

/**
 * @brief 销毁线程池，释放所有资源
 */
void thread_pool_destroy(ThreadPool* pool);

/**
 * @brief 单动量下并行批量模拟（多核压榨算力）
 * @param k 当前动量
 * @return 该动量下统计概率结果
 */
SingleKResult batch_run_single_k_parallel(double k);

#endif