#ifndef TRIPLET_H
#define TRIPLET_H

class Triplet{

public:
    int id;
    double value1, value2;

    Triplet(){}

    Triplet(int id, double value1, double value2): id(id), value1(value1), value2(value2) {}

    inline bool operator==(const Triplet & other) const {
        return (value1 == other.value1) && (value2 == other.value2);
    };

    inline bool operator!=(const Triplet & other) const {
        return (value1 != other.value1) || (value2 != other.value2);
    };

    inline bool operator<(const Triplet & other) const {
        if(value1 < other.value1){
            return true;
        } else if (value1 == other.value1){
            return value2 < other.value2;
        } else {
            return false;
        }
    }

    inline bool operator<=(const Triplet & other) const {
        if(value1 < other.value1){
            return true;
        } else if (value1 == other.value1){
            return value2 <= other.value2;
        } else {
            return false;
        }
    }

    inline bool operator>(const Triplet & other) const {
        if(value1 > other.value1){
            return true;
        } else if (value1 == other.value1){
            return value2 > other.value2;
        } else {
            return false;
        }
    }

    inline bool operator>=(const Triplet & other) const {
        if(value1 > other.value1){
            return true;
        } else if (value1 == other.value1){
            return value2 >= other.value2;
        } else {
            return false;
        }
    }
};

#endif //TRIPLET_H
