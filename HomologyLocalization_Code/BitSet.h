#ifndef BIT_SET_H
#define BIT_SET_H

#include <cstring>
#include <ostream>
#include <stdexcept>

class BitSet {
public:
	// Empty bit set
	BitSet() = default;

	// Allocating bit array. The BitSet is initialized with 0s
	BitSet(const int size) 
	{
		if (size > 0)
		{
			mAllocatedSize = ((size - 1)>> 3) + 1;
			mBitsize = size;
			mBits = (unsigned char *)calloc(mAllocatedSize, sizeof(unsigned char));
			if (mBits == nullptr) throw std::bad_alloc();
		}
	}

	~BitSet() 
	{
		if (mBits != nullptr) {
			free(mBits);
			mBits = nullptr;
		}
	}

	// Copy constructor
	BitSet(const BitSet & rhs) : mAllocatedSize(rhs.mAllocatedSize), mBitsize(rhs.mBitsize) 
	{
		if (mAllocatedSize != 0) 
		{
			if (mBits != nullptr)  // if not null, first release it
			{
				free(mBits);
				mBits = nullptr;
			}

			mBits = (unsigned char *)malloc(mAllocatedSize);
			memcpy(mBits, rhs.mBits, mAllocatedSize);
		}
	}

	// Copy constructor: copy the first 'bitCount' bits of rhs
	BitSet(const BitSet & rhs, int bitCount)
	{
		if (bitCount > 0)
		{
			if (bitCount > rhs.mBitsize)
				bitCount = rhs.mBitsize;

			mAllocatedSize = ((bitCount - 1) >> 3) + 1;
			mBitsize = bitCount;

			if (mBits != nullptr)  // if not null, first release it
			{
				free(mBits);
				mBits = nullptr;
			}

			mBits = (unsigned char *)malloc(mAllocatedSize);
			memcpy(mBits, rhs.mBits, mAllocatedSize);

			int numPackingBits =  mAllocatedSize * 8 - bitCount;
			mBits[mAllocatedSize - 1] <<= numPackingBits;
			mBits[mAllocatedSize - 1] >>= numPackingBits;
		}
	}

	BitSet & operator = (const BitSet & rhs) 
	{
		if (this != &rhs) {
			mAllocatedSize = rhs.mAllocatedSize;
			mBitsize = rhs.mBitsize;

			if (mBits != nullptr)  // if not null, first release it
			{
				free(mBits);
				mBits = nullptr;
			}

			mBits = (unsigned char *)malloc(mAllocatedSize);
			memcpy(mBits, rhs.mBits, mAllocatedSize);
		}
		return *this;
	}

	BitSet(BitSet && rhs) 
	{
		if (this == &rhs) { return; }
		if (mBits != nullptr)
		{
			free(mBits);
			mBits = nullptr;
		}

		mAllocatedSize = rhs.mAllocatedSize;
		mBits = rhs.mBits;
		mBitsize = rhs.mBitsize;
		rhs.mAllocatedSize = 0;
		rhs.mBits = nullptr;
		rhs.mBitsize = 0;
	}

	BitSet & operator = (BitSet && rhs) 
	{
		if (this == &rhs) { return  *this; }
		if (mBits != nullptr)  // if not null, first release it
		{
			free(mBits);
			mBits = nullptr;
		}

		mAllocatedSize = rhs.mAllocatedSize;
		mBits = rhs.mBits;
		mBitsize = rhs.mBitsize;
		rhs.mAllocatedSize = 0;
		rhs.mBits = nullptr;
		rhs.mBitsize = 0;
		return *this;
	}

	bool operator == (const BitSet & rhs) const
	{
		if (mBitsize != rhs.mBitsize)
			return false;

		if (memcmp(mBits, rhs.mBits, mAllocatedSize) != 0)
			return false;
		else
			return true;
	}

	// The following is about '<' and '>' operator
	bool operator < (const BitSet & rhs) const
	{
		if (mAllocatedSize < rhs.mAllocatedSize)
			return true;
		else if (mAllocatedSize > rhs.mAllocatedSize)
			return false;
		else if (memcmp(mBits, rhs.mBits, mAllocatedSize) < 0)
			return true;
		else 
			return false;
	}

	bool operator > (const BitSet & rhs) const
	{
		return rhs < *this;
	}

	// The following is about bit operator overriding: &, |, ^
	BitSet operator | (const BitSet & rhs)
	{
		if (mBitsize != rhs.mBitsize)
			throw std::invalid_argument(" The bits lengths of two operands are not the same! ");

		BitSet temp = *this;
		for (int i = 0; i < rhs.mAllocatedSize; i++)
			temp.mBits[i] |= rhs.mBits[i];

		return temp;
	}

	BitSet & operator |= (const BitSet & rhs)
	{
		if (mBitsize != rhs.mBitsize)
			throw std::invalid_argument(" The bits lengths of two operands are not the same! ");

		for (int i = 0; i < rhs.mAllocatedSize; i++)
			mBits[i] |= rhs.mBits[i];

		return *this;
	}

	BitSet operator & (const BitSet & rhs)
	{
		if (mBitsize != rhs.mBitsize)
			throw std::invalid_argument(" The bits lengths of two operands are not the same! ");

		BitSet temp = *this;
		for (int i = 0; i < rhs.mAllocatedSize; i++)
			temp.mBits[i] &= rhs.mBits[i];

		return temp;
	}

	BitSet & operator &= (const BitSet & rhs)
	{
		if (mBitsize != rhs.mBitsize)
			throw std::invalid_argument(" The bits lengths of two operands are not the same! ");

		for (int i = 0; i < rhs.mAllocatedSize; i++)
			mBits[i] &= rhs.mBits[i];

		return *this;
	}

	BitSet operator ^ (const BitSet & rhs)
	{
		if (mBitsize != rhs.mBitsize)
			throw std::invalid_argument(" The bits lengths of two operands are not the same! ");

		BitSet temp = *this;
		for (int i = 0; i < rhs.mAllocatedSize; i++)
			temp.mBits[i] ^= rhs.mBits[i];

		return temp;
	}

	BitSet & operator ^= (const BitSet & rhs)
	{
		if (mBitsize != rhs.mBitsize)
			throw std::invalid_argument(" The bits lengths of two operands are not the same! ");

		for (int i = 0; i < rhs.mAllocatedSize; i++)
			mBits[i] ^= rhs.mBits[i];

		return *this;
	}

	// Flip the bit at given position
	void flip(const int position)
	{
		if (position > mBitsize - 1 || position < 0)
			throw std::invalid_argument(" The positition to be set is beyond the scope! ");
		mBits[position >> 3] ^= (1 << (position & 0x7));
	}

	// Set the bit at given position to be 1
	void set(const int position) 
	{
		if (position > mBitsize - 1 || position < 0)
			throw std::invalid_argument(" The positition to be set is beyond the scope! ");
		mBits[position >> 3] |= (1 << (position & 0x7));
	}

	// Reset the bit at given position to be 0
	void reset(const int position)
	{
		if (position > mBitsize - 1 || position < 0)
			throw std::invalid_argument(" The positition to be set is beyond the scope! ");
		mBits[position >> 3] &= ~(1 << (position & 0x7));
	}

	// Check if the bit at given position is 0 or 1
	bool checkBit(const int position) const
	{
		if (position > mBitsize - 1 || position < 0)
			throw std::invalid_argument(" The positition to be set is beyond the scope! ");

		unsigned char temp = mBits[position >> 3];
		temp = temp << (7 - position & 0x7);
		temp = temp >> 7;

		return bool(temp);
	}

	// Set the bit at given position to a certain value
	void set(const size_t position, bool bitState)
	{
		if (bitState == true)
			set(position);
		else
			reset(position);
	}

	// Reset all the bits to be 0
	void reset()
	{
		memset(mBits, 0, sizeof(unsigned char) * mAllocatedSize);
	}

	// Return number of bits
	int getBitSize() const 
	{
		return mBitsize;
	}

	// Out Stream
	friend std::ostream& operator<<(std::ostream& out, const BitSet& bit_set) 
	{
		for (int i = bit_set.mAllocatedSize - 1; i >= 0; --i) {
			int byteLen = 8;
			unsigned char temp = '0';
			while (byteLen != 0)
			{
				temp = bit_set.mBits[i];
				temp = temp << (8 - byteLen);
				temp = temp >> 7;
				out << (int)temp;

				--byteLen;
			}
			out << " ";
		}
		return out;
	}

private:
	unsigned char* mBits = nullptr; // The array storing the bits
	
	int mAllocatedSize = 0; // Number of bytes allocated

	int mBitsize = 0; // Number of bits
};


#endif // !BIT_SET_H

