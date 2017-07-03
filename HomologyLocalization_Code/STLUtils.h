#ifndef STLUTILS_INCLUDED
#define STLUTILS_INCLUDED

template<typename ListT>
void mysort(ListT &l)
{
	std::sort(l.begin(), l.end());
}

template<typename ElemT>
void mysort(list<ElemT> &l)
{
	l.sort();
}

template<typename ListT>
void myclear(ListT &l)
{
	ListT().swap(l);
}

template<typename ElemT>
void myclear(list<ElemT> &l)
{
	l.clear();
}

template <class Container>
struct Counter
 : public std::iterator <std::output_iterator_tag,
                         void, void, void, void>
{ 
	size_t &cnt;

    Counter(size_t &x) : cnt(x) {}	
	Counter(const Counter&) = default;

	Counter& operator=(const Counter &) { return *this; }; // This line of code is necessary to pass the C++ 11 compiler
 
	template<typename t>
    Counter& operator=(t)
	{        
        return *this;
    }
    
    Counter& operator* () 
	{
        return *this;
    }
    
    Counter& operator++(int) 
	{
		++cnt;
        return *this;
    }    

	Counter& operator++() 
	{
		++cnt;
        return *this;
    }    
};


// We avoid excessive allocations by calculating the size of the resulting list.
// Then we resize the result and populate it with the actual values.
template<typename ListT>
ListT list_sym_diff(const ListT &sa, const ListT &sb){
	//assume inputs are both sorted increasingly	
	size_t count = 0;
	Counter<ListT> counter(count);
	set_symmetric_difference(sa.begin(), sa.end(), sb.begin(), sb.end(), counter);
	ListT out;
	
	out.reserve(count);
	set_symmetric_difference(sa.begin(), sa.end(), sb.begin(), sb.end(), back_inserter(out));	

	return out;
}
template<typename ListT>
ListT list_union(const ListT &sa, const ListT &sb){
	//assume inputs are both sorted increasingly	
	size_t count = 0;
	Counter<ListT> counter(count);
	set_union(sa.begin(), sa.end(), sb.begin(), sb.end(), counter);
	ListT out;
	
	out.reserve(count);
	set_union(sa.begin(), sa.end(), sb.begin(), sb.end(), back_inserter(out));	

	return out;
}
#endif
