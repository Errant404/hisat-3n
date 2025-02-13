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

#include "outq.h"

// beginRead：记录将要写入的读 id，并在 reorder 模式下扩展内部数组
void OutputQueue::beginRead(TReadId rdid, size_t threadId) {
	ThreadSafe t(&mutex_m, threadSafe_);
	nstarted_++;
	if(reorder_) {
		assert_geq(rdid, cur_);
		assert_eq(lines_.size(), finished_.size());
		assert_eq(lines_.size(), started_.size());
		if(rdid - cur_ >= lines_.size()) {
			// 扩展内部数组以确保能够存储当前读 id 的记录
			size_t oldsz = lines_.size();
			lines_.resize(rdid - cur_ + 1);
			started_.resize(rdid - cur_ + 1);
			finished_.resize(rdid - cur_ + 1);
			for(size_t i = oldsz; i < lines_.size(); i++) {
				started_[i] = finished_[i] = false;
			}
		}
		started_[rdid - cur_] = true;
		finished_[rdid - cur_] = false;
	}
}

// finishRead：完成当前读记录的写入
// 在 reorder 模式下，将记录存入内部数组；否则推入异步队列。
// 然后通过条件变量通知后台输出线程刷新数据。
void OutputQueue::finishRead(const BTString& rec, TReadId rdid, size_t threadId) {
	if(reorder_) {
		ThreadSafe t(&mutex_m, threadSafe_);
		assert_geq(rdid, cur_);
		assert_eq(lines_.size(), finished_.size());
		assert_eq(lines_.size(), started_.size());
		assert_lt(rdid - cur_, lines_.size());
		assert(started_[rdid - cur_]);
		assert(!finished_[rdid - cur_]);
		lines_[rdid - cur_] = rec;
		nfinished_++;
		finished_[rdid - cur_] = true;
	} else {
		// 非排序模式：将记录推入异步队列，不直接写出
		{
			std::unique_lock<std::mutex> lock(asyncMutex_);
			asyncQueue_.push(rec);
		}
		nfinished_++;
	}
	cv_.notify_one();
}

// flush：将缓冲中的记录写出到 obuf_
// 对于非排序模式，从 asyncQueue_ 中取出所有记录；
// 对于排序模式，依次写出连续的完成记录（达到阈值或强制刷新时）
void OutputQueue::flush(bool force, bool getLock) {
	if (!reorder_) {
		// 非排序模式：刷新异步队列
		std::unique_lock<std::mutex> lock(asyncMutex_);
		while(!asyncQueue_.empty()) {
			BTString rec = asyncQueue_.front();
			asyncQueue_.pop();
			obuf_.writeString(rec);
			nflushed_++;
		}
	} else {
		// 排序模式：获取 mutex_m 后检查连续完成的记录
		ThreadSafe t(&mutex_m, getLock && threadSafe_);
		size_t nflush = 0;
		while(nflush < finished_.size() && finished_[nflush]) {
			nflush++;
		}
		if(force || nflush >= NFLUSH_THRESH) {
			for(size_t i = 0; i < nflush; i++) {
				obuf_.writeString(lines_[i]);
			}
			// 擦除已写出的记录，更新当前读 id
			lines_.erase(0, nflush);
			started_.erase(0, nflush);
			finished_.erase(0, nflush);
			cur_ += nflush;
			nflushed_ += nflush;
		}
	}
}

// asyncOutputThread：后台输出线程
// 线程循环等待条件变量通知或超时，然后调用 flush() 将缓冲数据写出。
// 当 shutdown_ 标志为 true 时退出循环，并进行最后一次刷新。
void OutputQueue::asyncOutputThread() {
	while (true) {
		{
			std::unique_lock<std::mutex> lock(cv_mutex_);
			cv_.wait_for(lock, std::chrono::milliseconds(50), [this](){ return shutdown_; });
			if(shutdown_)
				break;
		}
		flush(true, true);
	}
	// 退出前进行最后一次刷新
	flush(true, true);
} 

#ifdef OUTQ_MAIN

#include <iostream>

using namespace std;

int main(void) {
	cerr << "Case 1 (one thread) ... ";
	{
		OutFileBuf ofb;
		OutputQueue oq(ofb, false);
		assert_eq(0, oq.numFlushed());
		assert_eq(0, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(1);
		assert_eq(0, oq.numFlushed());
		assert_eq(1, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(3);
		assert_eq(0, oq.numFlushed());
		assert_eq(2, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(2);
		assert_eq(0, oq.numFlushed());
		assert_eq(3, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(3, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(0);
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.finishRead(0);
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.flush(true);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.finishRead(2);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(2, oq.numFinished());
		oq.flush(true);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(2, oq.numFinished());
		oq.finishRead(1);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(3, oq.numFinished());
		oq.flush(true);
		assert_eq(3, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(3, oq.numFinished());
	}
	cerr << "PASSED" << endl;

	cerr << "Case 2 (one thread) ... ";
	{
		OutFileBuf ofb;
		OutputQueue oq(ofb, false);
		BTString& buf1 = oq.beginRead(0);
		BTString& buf2 = oq.beginRead(1);
		BTString& buf3 = oq.beginRead(2);
		BTString& buf4 = oq.beginRead(3);
		BTString& buf5 = oq.beginRead(4);
		assert_eq(5, oq.numStarted());
		assert_eq(0, oq.numFinished());
		buf1.install("A\n");
		buf2.install("B\n");
		buf3.install("C\n");
		buf4.install("D\n");
		buf5.install("E\n");
		oq.finishRead(4);
		oq.finishRead(1);
		oq.finishRead(0);
		oq.finishRead(2);
		oq.finishRead(3);
		oq.flush(true);
		assert_eq(5, oq.numFlushed());
		assert_eq(5, oq.numStarted());
		assert_eq(5, oq.numFinished());
		ofb.flush();
	}
	cerr << "PASSED" << endl;
	return 0;
}

#endif /*def ALN_SINK_MAIN*/
