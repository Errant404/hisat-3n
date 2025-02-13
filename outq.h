/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OUTQ_H_
#define OUTQ_H_

#include "assert_helpers.h"
#include "ds.h"
#include "sstring.h"
#include "read.h"
#include "threading.h"
#include "mem_ids.h"

#include <thread>
#include <condition_variable>
#include <queue>
#include <mutex>
#include <chrono>

class OutputQueue {

    static const size_t NFLUSH_THRESH = 8;

public:
    // 构造函数接口保持不变。参数含义：
    // obuf      ：用于写出结果的缓冲区
    // reorder   ：是否需要对输出进行重新排序（按读 id 顺序）
    // nthreads  ：使用线程数（当 nthreads>1 时要求 threadSafe 为 true）
    // threadSafe：是否启用线程安全机制
    // rdid      ：起始读 id
    OutputQueue(OutFileBuf& obuf, bool reorder, size_t nthreads, bool threadSafe, TReadId rdid = 0) :
        obuf_(obuf),
        cur_(rdid),
        nstarted_(0),
        nfinished_(0),
        nflushed_(0),
        lines_(RES_CAT),
        started_(RES_CAT),
        finished_(RES_CAT),
        reorder_(reorder),
        threadSafe_(threadSafe),
        mutex_m(),
        shutdown_(false)
    {
        assert(nthreads <= 1 || threadSafe);
        // 启动后台异步输出线程
        output_thread_ = std::thread(&OutputQueue::asyncOutputThread, this);
    }

    // 析构函数：退出前通知异步线程结束并等待线程退出
    ~OutputQueue() {
        {
            std::unique_lock<std::mutex> lock(cv_mutex_);
            shutdown_ = true;
        }
        cv_.notify_all();
        if (output_thread_.joinable())
            output_thread_.join();
    }

    // 调用 beginRead 告知队列将要写入读 id 为 rdid 的记录
    void beginRead(TReadId rdid, size_t threadId);

    // 调用 finishRead 完成输出记录的写入。注意：在异步模式下，该方法仅将记录推入队列，
    // 后台线程负责实际写出，减少锁竞争。
    void finishRead(const BTString& rec, TReadId rdid, size_t threadId);

    // 返回当前缓冲的记录数（对于 reorder 模式返回内部数组大小，
    // 对于非排序模式返回异步队列大小）
    size_t size() const {
        return reorder_ ? lines_.size() : asyncQueue_.size();
    }
    
    // 返回已刷新（写出）的记录数
    TReadId numFlushed() const { return nflushed_; }
    // 返回已开始写入的记录数
    TReadId numStarted() const { return nstarted_; }
    // 返回已完成写入的记录数
    TReadId numFinished() const { return nfinished_; }

    // flush 用于将缓冲中的记录写出到 obuf_。force 为 true 时强制刷新。
    // getLock 指示是否获取内部锁。
    void flush(bool force = false, bool getLock = true);

protected:
    // 异步输出线程执行函数。该线程不断等待通知或超时，然后调用 flush() 进行刷新。
    void asyncOutputThread();

    OutFileBuf&     obuf_;
    TReadId         cur_;
    TReadId         nstarted_;
    TReadId         nfinished_;
    TReadId         nflushed_;
    EList<BTString> lines_;    // 用于存储排序模式下的输出记录
    EList<bool>     started_;  // 标记每个记录是否已开始写入
    EList<bool>     finished_; // 标记每个记录是否已完成写入
    bool            reorder_;  // 是否需要按照读 id 顺序输出
    bool            threadSafe_;
    MUTEX_T         mutex_m;   // 保护 reorder_ 模式下数据的锁

    // 针对非 reorder 模式，采用一个简单的异步队列存储输出记录
    std::queue<BTString> asyncQueue_;
    std::mutex asyncMutex_;

    // 异步输出线程及其条件变量，用于唤醒刷新操作
    std::thread output_thread_;
    std::condition_variable cv_;
    std::mutex cv_mutex_;
    bool shutdown_;
};


// 辅助类：构造时调用 beginRead，析构时调用 finishRead
class OutputQueueMark {
public:
    OutputQueueMark(OutputQueue& q, const BTString& rec, TReadId rdid, size_t threadId) :
        q_(q),
        rec_(rec),
        rdid_(rdid),
        threadId_(threadId)
    {
        q_.beginRead(rdid, threadId);
    }
    
    ~OutputQueueMark() {
        q_.finishRead(rec_, rdid_, threadId_);
    }
    
protected:
    OutputQueue& q_;
    const BTString& rec_;
    TReadId rdid_;
    size_t threadId_;
};

#endif
