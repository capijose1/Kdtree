#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

//bounded
using namespace std;
template <size_t N, typename ElemType>
class KDTree {
public:
    typedef std::pair<Point<N>, ElemType> value_type;
    KDTree();
    ~KDTree();
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    size_t dimension() const;
    size_t size() const;
    bool empty() const;
    bool contains(const Point<N>& pt) const;
    void insert(const Point<N>& pt, const ElemType& value = ElemType());
    ElemType& operator[](const Point<N>& pt);
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    ElemType knn_value(const Point<N>& key, size_t k) const;
    vector<ElemType> knn_query(const Point<N>& key, size_t k) const {
        std::vector<ElemType> values;
        return values;
    }

private:
    size_t dimension_;
    size_t size_;
    vector<value_type> data;
    struct Nodo {
        Nodo() = default;
        Point<N> point;
        int nivel;
        ElemType valor;
        Nodo* izquierda;
        Nodo* derecha;
        Nodo(const Point<N>& ptr, int level, const ElemType& value = ElemType()) {
            point = ptr;
            izquierda = derecha = NULL;
            nivel = level;
            valor = value;

        }

    };
    Nodo* root_;
    Nodo* copiar(Nodo* root) {
        if (root == NULL) return NULL;
        Nodo* nod = new Nodo(*root);
        nod->izquierda = copiar(root->izquierda);
        nod->derecha = copiar(root->derecha);
        return nod;
    }
    void limpiar(Nodo* root) {
        if (root == NULL) return;
        limpiar(root->izquierda);
        limpiar(root->derecha);
        delete root;
    }
    Nodo* buscar(typename KDTree<N, ElemType>::Nodo* nodo, const Point<N>& pt) const {
        if (nodo == NULL || nodo->point == pt) return nodo;
        int nivel = nodo->nivel;
        if (pt[nivel % N] < nodo->point[nivel % N]) {
            if (nodo->izquierda) {
                return buscar(nodo->izquierda, pt);
            }
            return nodo;
        }
        else {
            if (nodo->derecha) {
                return buscar(nodo->derecha, pt);
            }
            return nodo;
        }
    }

};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    size_ = 0;
    root_ = NULL;
    dimension_ = N;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    limpiar(root_);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
    root_ = copiar(rhs.root_);
    size_ = rhs.size_;
    dimension_ = rhs.dimension_;
}
template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) { //copiar nuevo arbol
    if (this != &rhs) { // el this es el arbol a comparar 
        root_ = copiar(rhs.root_);
        size_ = rhs.size_;
        dimension_ = rhs.dimension_;
    }
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    return size_ == 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
    return  buscar(root_, pt)->point == pt;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    Nodo* nodo = buscar(root_, pt);
    if (nodo == NULL) {
        root_ = new Nodo(pt, 0, value);
        size_ = 1;
    }
    else {
        if ((nodo)->point == pt) {
            (nodo)->valor = value;
        }
        else {
            int nivel = (nodo)->nivel;
            Nodo* newNodo = new Nodo(pt, nivel + 1, value);
            if (pt[nivel % N] < (nodo)->point[nivel % N]) {
                (nodo)->izquierda = newNodo;
            }
            else {
                (nodo)->derecha = newNodo;
            }
            size_++;
        }

    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
    Nodo* node = buscar(root_, pt);
    if (node != NULL && node->point == pt) {
        return node->valor;
    }
    else { 
        insert(pt);
        if (node == NULL) {
            return root_->valor;
        }
        else {
            if (node->izquierda != NULL && node->izquierda->point == pt) {
                return node->izquierda->valor;
            }
            else {
                return node->derecha->valor;
            }
        }
    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
    return const_cast<ElemType&>(const_cast<const KDTree&>(*this).at(pt));
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    Nodo* n = buscar(root_, pt);
    if (n == NULL || n->point != pt) {
        throw  out_of_range("Fuera de rango");
    }
    return (n)->valor;
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N>& key, size_t k) const {
    ElemType new_element;
    return new_element;
}



// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // S

//ayuda de https://github.com/robertoBenavides/kdtree/blob/main/src/KDTree.hpp , https://rosettacode.org/wiki/K-d_tree