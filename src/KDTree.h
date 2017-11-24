/**
 * File: KDTree.h
 * Author: (your name here)
 * ------------------------
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 */

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <queue>
#include <map>
using namespace std;

template <size_t N, typename ElemType>
class KDTree {
public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();
    
    // Destructor: ~KDTree()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTree.
    ~KDTree();
    
    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Deep-copies the contents of another KDTree into this one.
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;
    
    // size_t size() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree and whether the tree is
    // empty.
    size_t size() const;
    bool empty() const;
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>& pt) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);
    
    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    
    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;

private:
    struct node {
        Point<N> point;
        ElemType elem;
        node* left;
        node* right;
    };
    node* start_node;
    size_t size_count;
    void deleteNodes(node *&n);
    node* findNode(const Point<N> & pt, int & found_state) const;
    void addLeaf( node* & n, const int & found_state,const Point<N> & pt, const ElemType & value);

    void traverse(const node *n, size_t level, const Point<N> &pt,
                  BoundedPQueue<ElemType> & pq_dist) const;
    void copyNodes(node* & copied, const node *other);
    void copy(const KDTree& rhs);
    void clear();
    void initialise();

};

/** KDTree class implementation details */

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree():start_node(new node),size_count(0) {

}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::deleteNodes(node * & n) {
    if(size_count>0){
        if(n->left==NULL && n->right==NULL){
            delete n;
        }else {
            if(n->left!=NULL){
                deleteNodes(n->left);
            }
            if(n->right!=NULL){
                deleteNodes(n->right);
            }
        }
    }

}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::clear(){
    deleteNodes(start_node);
    size_count = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    clear();
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>:: size() const{
    return size_count;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const{
    return size_count == 0;

}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>:: node * KDTree<N, ElemType>::
findNode(const Point<N> & pt, int &found_state) const{
    node * curr_node=start_node;
    Point<N> curr_point;
    size_t rem;
    found_state=0;

    if(size_count==0){
        return start_node;
    }

    for (size_t i = 0; i < size_count; i++){
        if(curr_node!=NULL){
            rem= i % N;
            curr_point=curr_node->point;
            if(pt == curr_point){
                found_state=1;
                break;
            }
            if(pt[rem] <= curr_point[rem]){
                if(curr_node->left==NULL){
                    found_state=-1;
                    break;
                }else{
                    curr_node=curr_node->left;
                }

            }else{
                if(curr_node->right==NULL){
                    found_state=2;
                    break;
                }else{
                    curr_node=curr_node->right;
                }
            }
        }else{
            break;
        }
    }
    return curr_node;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>:: addLeaf( node* & n, const int & found_state,const Point<N> & pt, const ElemType & value){
    node* new_node = new node;
    new_node->point=pt;
    new_node->elem=value;
    new_node->left=NULL;
    new_node->right=NULL;
    switch (found_state) {
    case 0:
        *n= *new_node;
        break;
    case -1:
        n->left = new_node;
        break;
    case 2:
        n->right = new_node;
    default:
        break;
    }
    size_count++;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert (const Point<N>& pt, const ElemType& value){
    int found_state;
    node* found_node = findNode(pt, found_state);
    if(found_state==1){
        found_node->elem=value;
    }else{
        addLeaf(found_node,found_state,pt,value);
    }

}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>:: contains(const Point<N>& pt) const{
    int found_state;
    findNode(pt,found_state);
    return found_state==1;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>:: operator[] (const Point<N>& pt){
    int found_state;
    node* found_node= findNode(pt,found_state);

    if(found_state!=1){
        ElemType value = numeric_limits<ElemType>::max();
        addLeaf(found_node,found_state,pt,value);
        found_node=findNode(pt,found_state);
    }
    return (found_node->elem);
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>:: at(const Point<N>& pt){
    int found_state;
    node* found_node= findNode(pt,found_state);

    if(found_state==1){
        return (found_node->elem);
    }else{
        throw out_of_range("Invalid point.");
    }
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>:: at(const Point<N>& pt) const{
    int found_state;
    node* found_node= findNode(pt,found_state);

    if(found_state==1){
        return (found_node->elem);
    }else{
        throw out_of_range("Invalid point.");
    }
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::traverse(const node* n, size_t level, const Point<N> &pt,
                                   BoundedPQueue<ElemType> &pq_dist) const{
    Point<N> curr_point;
    size_t rem;
    double distance;

    if(n!=NULL){
        rem= level % N;
        level++;
        curr_point=n->point;
        distance=Distance(curr_point,pt);
        pq_dist.enqueue(n->elem,distance);

        if(pt[rem] <= curr_point[rem]){
            if(n->left!=NULL){
                traverse(n->left,level,pt,pq_dist);
                if(pq_dist.size()<pq_dist.maxSize() || fabs(pt[rem]-curr_point[rem]) < pq_dist.best()){
                    traverse(n->right,level,pt,pq_dist);
                }
            }
        }else{
            if(n->right!=NULL){
                traverse(n->right,level,pt,pq_dist);
                if(pq_dist.size()<pq_dist.maxSize() || abs(pt[rem]-curr_point[rem]) < pq_dist.best()){
                    traverse(n->left,level,pt,pq_dist);
                }
            }
        }
    }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const{

    BoundedPQueue<ElemType> pq_dist(k);
    if(size_count==0){
        throw out_of_range("The k-D tree is empty.");
    }else{
        traverse(start_node,0,key,pq_dist);
    }

    ElemType elem_best;
    ElemType elem_curr;
    size_t freq_max=0;
    size_t freq;
    map<ElemType,size_t> map_freq;

    while(!pq_dist.empty()){
        elem_curr=pq_dist.dequeueMin();
        if(map_freq.count(elem_curr)==0){
            map_freq.insert(make_pair(elem_curr,1));
        }else{
            map_freq[elem_curr]++;
        }
        freq=map_freq[elem_curr];
        if(freq > freq_max){
            freq_max=freq;
            elem_best=elem_curr;
        }
    }
    return elem_best;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::copyNodes(node *&copied, const node *other){

    copied->point = other->point;
    copied->elem=other->elem;
    if(other->left!=NULL){
        copied->left=new node;
        copyNodes(copied->left,other->left);
    }else {
        copied->left=NULL;
    }

    if(other->right!=NULL){
        copied->right=new node;
        copyNodes(copied->right,other->right);
    }else{
        copied->right=NULL;
    }
    return;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::copy(const KDTree& rhs){


    if (rhs.size()!=0){
        copyNodes(start_node,rhs.start_node);
        size_count=rhs.size();
    }

}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs):start_node(new node),size_count(0){
    if(&rhs!=this){

        copy(rhs);
    }
}

template <size_t N, typename ElemType>
KDTree<N,ElemType>& KDTree<N, ElemType>:: operator=(const KDTree& rhs){
    if(&rhs!=this){

        clear();
        start_node= new node;
        copy(rhs);
    }
    return *this;
}

#endif // KDTREE_INCLUDED
