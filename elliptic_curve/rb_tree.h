#ifndef RB_TREE_H 
#define RB_TREE_H 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "ec_point.h"

typedef ec_point T;
#define compLT(a,b) ec_compare_bs(a,b)<0
#define compEQ(a,b) ec_compare_bs(a,b)==0

typedef enum { BLACK, RED } nodeColor;

typedef struct Node_ {
    struct Node_ *left;       
    struct Node_ *right;      
    struct Node_ *parent;     
    nodeColor color;
    T data;
} Node;

class rb_tree {

private:
Node * root;
Node sentinel;
void rotateLeft(Node *x);
void rotateRight(Node *x); 
void insertFixup(Node *x); 

void deleteFixup(Node *x);
void deleteNode(Node *z);

public:
Node *insertNode(T data);
Node *findNode(T data);
rb_tree();
};

#define NIL &sentinel

#endif
