//
// Created by pentagon03 on 8/9/24.
//

#ifndef ALIGNASM_PRIORITY_QUEUE_VECTOR_HPP
#define ALIGNASM_PRIORITY_QUEUE_VECTOR_HPP
#include<queue>
#include<vector>
#include<utility>
#include<algorithm>

template<typename Tp, typename Sequence = std::vector<Tp>, typename Compare  = std::less<typename Sequence::value_type> >
class PQVec : public std::priority_queue<Tp, Sequence, Compare> {
public:
    using Vec = std::vector<Tp>;
    const Vec& getVector(){return this->c;}
    Vec getSortedVector() {
        Vec r(this->c.begin(),this->c.end());
        // c is a heap
        std::sort_heap(r.begin(), r.end(), this->comp);
        // priority-queue order
        std::reverse(r.begin(), r.end());
        return std::move(r);
    }
};
#endif //ALIGNASM_PRIORITY_QUEUE_VECTOR_HPP
