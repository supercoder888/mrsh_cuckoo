/* 
 * File:   util.h
 * Author: mustafakarabat
 *
 * Created on 5. Juni 2012, 14:38
 */
#ifndef HELPER_H
#define	HELPER_H
#include <stdio.h>
#include <stdint.h>
#include "../header/config.h"


FILE    		*getFileHandle(char *filename);
unsigned int	find_file_size(FILE *fh);



#endif	/* HELPER_H */

