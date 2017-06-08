#ifndef VECTORELEMENT_H
#define VECTORELEMENT_H

class VectorElement {
public:
	int id;
	double data;

	inline VectorElement() { }

	inline VectorElement(int _id, double _data) {
		id = _id;
		data = _data;
	}

	inline bool operator==(const VectorElement & other) const {
		if (this == &other)
			return true;
		return data == other.data && id == other.id;
	};

	inline bool operator!=(const VectorElement & other) const {
		if (this == &other)
			return false;
		return data != other.data || id != other.id;
	};

	inline bool operator<(const VectorElement & other) const {
		return data < other.data;
	}

	inline bool operator<=(const VectorElement & other) const {
		return data <= other.data;
	}

	inline bool operator>(const VectorElement & other) const {
		return data > other.data;
	}

	inline bool operator>=(const VectorElement & other) const {
		return data >= other.data;
	}

};

#endif //VECTORELEMENT_H
