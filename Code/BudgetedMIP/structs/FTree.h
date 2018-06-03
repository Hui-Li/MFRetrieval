#ifndef FTREE_H
#define FTREE_H

#include "../util/Base.h"

// An implementation of Algorithm 3: selection tree described in TAOCP by Knuth [8]
template<typename T>
struct Ftree_t { // {{{
	size_t len ;     // real length of valid elements
	size_t elements; // 2^ceil(log2(len)) capacity
	std::vector<T> val;
	T *true_val;

	Ftree_t(size_t len_= 0) {resize(len_);}
	Ftree_t(const Ftree_t &other) { // {{{
		len = other.len;
		elements = other.elements;
		val = other.val;
		true_val = &val[elements];
	} // }}}
	Ftree_t& operator=(const Ftree_t &other) { // {{{
		len = other.len;
		elements = other.elements;
		val = other.val;
		true_val = &val[elements];
	} // }}}k
	size_t size() const {return len;}
	void resize(size_t len_) { // {{{
		len = len_;
		if(len == 0){
			val.clear(); elements = 0; true_val = NULL;
			return;
		} else {
			elements = 1;
			while (elements<len) elements <<= 1;
			val.clear(); val.resize(2*elements, T::default_value());
			true_val = &val[elements];
		}
	} // }}}
	void clear() {resize(0);}

	T &operator[](size_t idx) {return true_val[idx];}
	const T &operator[](size_t idx) const {return true_val[idx];}

	void init_dense() { // {{{
		for(size_t pos = elements-1; pos > 0; --pos)
			val[pos] = T::op(val[pos<<1], val[(pos<<1)+1]);
	} // }}}
	void set_value(size_t idx, const T& v) { // {{{
		val[idx+=elements] = v;
		while (idx>1) {
			val[idx>>1] = T::op(val[idx],val[idx^1]);
			idx >>=1;
		}
	} // }}}
	void print() { // {{{
		std::cout << "len " << len << "; elements " << elements << std::endl;
		size_t idx = 1;
		while (idx <= elements) {
			for(size_t pos = idx; pos < (idx << 1); pos++)
				std::cout << "(" << val[pos].idx<<"," <<val[pos].value<< ") ";
			std::cout << std::endl;
			idx <<= 1;
		}
	} // }}}
}; // }}}

#endif //FTREE_H
