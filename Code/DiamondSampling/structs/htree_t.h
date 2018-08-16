#ifndef HTREE_T_H
#define HTREE_T_H

#include "../util/Base.h"

template<typename T=double>
struct htree_t{ // {{{
	size_t size;     // real size of valid elements
	size_t elements; // 2^ceil(log2(size)) capacity
	std::vector<T> val;
	T *true_val;

	T& operator[](size_t idx) { assert(idx < elements); return true_val[idx]; }
	const T& operator[] (size_t idx) const { assert(idx < elements); return true_val[idx]; }
	void init_dense() { // {{{
		/*
		for(size_t pos = (elements+size)>>1; pos > 0; --pos)
			val[pos] = val[pos<<1] + val[(pos<<1)+1];
			*/
		for(size_t pos = elements-1; pos > 0; --pos)
			val[pos] = val[pos<<1] + val[(pos<<1)+1];
	} // }}}
	void update_parent(size_t idx, T delta) { // {{{
		idx = (idx+elements)>>1;
		while(idx) {
			val[idx] += delta;
			idx >>= 1;
		}
	} // }}}
	void set_value(size_t idx, T value) { // {{{
		value -= val[idx+=elements]; // delta update
		while(idx) {
			val[idx] += value;
			idx >>= 1;
		}
	} // }}}
	// urnd: uniformly random number between [0,1]
	size_t log_sample(double urnd) { // {{{
		//urnd *= val[1];
		size_t pos = 1;
		while(pos < elements) {
		//while(pos < size) {
			pos <<= 1;
			if(urnd > val[pos])
				urnd -= val[pos++];
			/*
			double tmp = urnd - val[pos];
			if(tmp >= 0) {
				urnd = tmp;
				pos++;
			}
			*/
			/*
			if(urnd < val[pos*2])
				pos = pos*2;
			else {
				urnd -= val[pos*2];
				pos = pos*2+1;
			}
			*/
		}
		return pos-elements;
	} // }}}
	size_t linear_sample(double urnd) { // {{{
		//urnd = urnd*val[1];
		size_t pos = elements;
		while(urnd > 0)
			urnd -= val[pos++];
		if(pos >= elements+size) pos = elements+size-1;
		return pos-elements;
	} // }}}
	double total_sum() { return val[1]; }
	double left_cumsum(size_t idx) { // {{{
		if(idx == elements) return val[1];
		size_t pos = elements+idx+1;
		double sum = 0;
		while(pos>1) {
			if(pos & 1)
				sum += val[pos^1];
			pos >>= 1;
		}
		return sum;
	} // }}}
	double right_cumsum(size_t idx) {return val[1] - left_cumsum(idx-1);}
	htree_t(size_t size=0) { resize(size); }
	htree_t(const htree_t& other) {
		size = other.size;
		elements = other.elements;
		val = other.val;
		true_val = &val[elements];
	}
	htree_t& operator=(const htree_t& other) {
		size = other.size;
		elements = other.elements;
		val = other.val;
		true_val = &val[elements];
	}

	void resize(size_t size_) { // {{{
		size = size_;
		if(size == 0) {
			val.clear(); elements = 0; true_val = NULL;
			return;
		} else {
			elements = 1;
			while(elements < size) elements <<=1;
			val.clear(); val.resize(2*elements, 0);
			true_val = &val[elements];
		}
		//Q.reserve(elements);
	} //}}}
	void clear() { for(auto &v: val) v = 0; }
}; // }}}

#endif //HTREE_T_H
