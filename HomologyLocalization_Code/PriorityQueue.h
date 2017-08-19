#ifndef _PRIORITY_QUEUE_H_
#define _PRIORITY_QUEUE_H_
#include <iostream>
#include <map>
#include <vector>
#include "Globals.h"
#include "BitSet.h"
using namespace std;


class priorityQueue
{
public:
	priorityQueue();
	~priorityQueue();
	void push(const cgNode &);
	cgNode top();
	void pop();
	void clear();
	bool isEmpty();
	bool updateNodeInQueue(const cgNode &);

private:
	vector<cgNode> dataArray;
	int mSize;
	void heapAdjustPop(int, int);
	void heapAdjustPush(int);
	map<std::pair<int, BitSet>, int> mMap;
};


priorityQueue::priorityQueue()
{
	mSize = 0; // the 0th element will not be used
	dataArray.push_back(cgNode()); // push a place holder
}

priorityQueue::~priorityQueue()
{
	dataArray.clear();
	mMap.clear();
}

void priorityQueue::push(const cgNode & newItem)
{
	if (mSize == (dataArray.capacity() - 1)) // no enough memory
	{
		dataArray.resize(dataArray.capacity() * 2);
	}
	++mSize;
	dataArray[mSize] = newItem;
	mMap[std::make_pair(newItem.vertex, newItem.sumAnnotation)] = mSize;

	heapAdjustPush(mSize);
}

cgNode priorityQueue::top()
{
	return dataArray[1];
}

void priorityQueue::pop()
{
	if (mSize == 0)
	{
		cout << "Cannot pop any more!" << endl;
		system("pause");
		exit(1);
	}
	cgNode tempFront = dataArray[1];
	dataArray[1] = dataArray[mSize--];
	heapAdjustPop(1, mSize);
	mMap.erase(std::make_pair(tempFront.vertex, tempFront.sumAnnotation));
}

void priorityQueue::clear()
{
	mSize = 0; // release the memory
	dataArray.clear();
	mMap.clear();

	dataArray.push_back(cgNode()); // push a place holder
}

void priorityQueue::heapAdjustPop(int start, int end) // the start element is to be sifted
{
	int s = start;
	cgNode rc = dataArray[start];
	for (int j = 2 * s; j <= end; j *= 2)
	{
		if (j < end && (dataArray[j] < dataArray[j + 1]))
			++j;
		if (!(rc < dataArray[j]))
			break;

		dataArray[s] = dataArray[j];
		mMap[std::make_pair(dataArray[j].vertex, dataArray[j].sumAnnotation)] = s;
		s = j;
	}
	dataArray[s] = rc;
	mMap[std::make_pair(rc.vertex, rc.sumAnnotation)] = s;
}

void priorityQueue::heapAdjustPush(int end)
{
	int parent = end / 2;
	int curr = end;
	cgNode temp_curr, temp_parent;
	while (parent != 0)
	{
		if (dataArray[parent] < dataArray[curr])
		{
			temp_parent = dataArray[parent];
			temp_curr = dataArray[curr];
			dataArray[parent] = temp_curr;
			dataArray[curr] = temp_parent;

			mMap[std::make_pair(temp_curr.vertex, temp_curr.sumAnnotation)] = parent;
			mMap[std::make_pair(temp_parent.vertex, temp_parent.sumAnnotation)] = curr;

			curr = parent;
			parent = parent / 2;
		}
		else
		{
			break;
		}
	}
}

bool priorityQueue::isEmpty()
{
	return mSize == 0;
}

bool priorityQueue::updateNodeInQueue(const cgNode & newNode)
{
	map<pair<int, BitSet>, int>::iterator it;

	it = mMap.find(std::make_pair(newNode.vertex, newNode.sumAnnotation));
	if (it == mMap.end())
	{
		return false;
	}
	else
	{
		int currGScore = dataArray[it->second].gScore;
		if (currGScore > newNode.gScore)
		{
			int hVal = dataArray[it->second].fScore - dataArray[it->second].gScore;
			dataArray[it->second].gScore = newNode.gScore;
			dataArray[it->second].fScore = hVal + newNode.gScore;
			dataArray[it->second].previous = newNode.previous;

			// remake the heap
			heapAdjustPush(it->second);
		}
		return true;
	}
}

#endif