#ifndef FASTHEAP_H
#define FASTHEAP_H

typedef double datatype;

/*returns the parent of a heap position*/
inline static int parent(int posel) {
    if (posel % 2) /*odd*/
        return posel / 2;
    else
        return (posel - 1) / 2;
}

/*enqueues element*/
/*heap[0] is the element with the smallest data*/
/*every element is greater than or equal to its parent*/
inline static void heap_enqueue(datatype el, int id, VectorElement *heap, int *num_elems) {
    datatype tmp;
    int tmpidx;
    int p;
    int posel;

    posel = *num_elems; //last position
    heap[(*num_elems)].data = el;
    heap[(*num_elems)++].id = id;

    while (posel > 0) {
        p = parent(posel);
        if (el < heap[p].data) {
            /* swap element with its parent */
            tmp = heap[p].data;
            heap[p].data = el;
            heap[posel].data = tmp;

            tmpidx = heap[p].id;
            heap[p].id = id;
            heap[posel].id = tmpidx;

            posel = parent(posel);
        }
        else break;
    }
}

/* moves down the root element */
/* used by dequeue (see below) */
inline static void heap_movedown(VectorElement *heap, int *num_elems) {
    datatype tmp;
    int tmpidx;
    int posel = 0; //root
    int swap;
    /*while posel is not a leaf and heap[posel].data > any of childen*/
    while (posel * 2 + 1 < *num_elems) /*there exists a left son*/
    {
        if (posel * 2 + 2 < *num_elems) /*there exists a right son*/
        {
            if (heap[posel * 2 + 1].data < heap[posel * 2 + 2].data)
                swap = posel * 2 + 1;
            else
                swap = posel * 2 + 2;
        }
        else
            swap = posel * 2 + 1;

        if (heap[posel].data > heap[swap].data) /*larger than smallest son*/
        {
            /*swap elements*/
            tmp = heap[swap].data;
            heap[swap].data = heap[posel].data;
            heap[posel].data = tmp;

            tmpidx = heap[swap].id;
            heap[swap].id = heap[posel].id;
            heap[posel].id = tmpidx;

            posel = swap;
        }
        else break;
    }
}

/* returns the root element, puts the last element as root and moves it down */
inline static void heap_dequeue(VectorElement *heap, int *num_elems) {
    if ((*num_elems) == 0) /* empty queue */
        return;

//    el->data = heap[0].data;
//    el->id = heap[0].id;
    heap[0].data = heap[(*num_elems) - 1].data;
    heap[0].id = heap[(*num_elems) - 1].id;
    (*num_elems)--;
    heap_movedown(heap, num_elems);
}

/* prints the elements in the heap (starting from the top element)*/
inline void heap_print(VectorElement *heap, int num_elems) {

    for (int i = 0; i < num_elems; i++)
        cout << heap[i].data << "," << heap[i].id << endl;

}

#endif //FASTHEAP_H
