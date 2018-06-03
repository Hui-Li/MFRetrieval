#ifndef COMPARATOR_H
#define COMPARATOR_H

#include "Base.h"

// descending comparator for arg_sort
template<typename T>
struct comparator {
	const T* value;
	comparator(const T* value=NULL): value(value){}
	bool operator()(size_t a, size_t b) const { return value[a] > value[b]; }
};

#endif //COMPARATOR_H
