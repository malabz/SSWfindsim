#ifndef __UTILS_H__
#define __UTILS_H__
#include "sswfindsim.h"

typedef struct _link_list
{
	struct _link_list* nxt;
	char* s;
} link_list;

link_list* new_node();
void free_node(link_list* this_node);

double countATGC(char* s);

void alignmentout(s_align* a, FILE *filename, int number);

#endif