#ifndef LUTRO_ROTOZOOM_H
#define LUTRO_ROTOZOOM_H

#include <stdint.h>
#include <math.h>

#include "painter.h"


#ifndef M_PI
#define M_PI	3.1415926535897932384626433832795
#endif

/* Set up for C function definitions, even when using C++ */
#ifdef __cplusplus
extern "C" {
#endif
	/* ---- Defines */

	/*!
	\brief Disable anti-aliasing (no smoothing).
	*/
#define SMOOTHING_OFF		0

	/*!
	\brief Enable anti-aliasing (smoothing).
	*/
#define SMOOTHING_ON		1

	/* ---- Function Prototypes */

	/* 

	Rotozoom functions

	*/

	bitmap_t * transform_rotozoom(const bitmap_t * src, double angle, double zoom, int smooth);

	bitmap_t * transform_rotozoom_xy(const bitmap_t * src, double angle, double zoomx, double zoomy, int smooth);


	void get_rotozoom_size(int width, int height, double angle, double zoom, int *dstwidth, int * dstheight);

	void get_rotozoom_size_xy(
		int width,
		int height,
		double angle,
		double zoomx,
		double zoomy, 
		int *dstwidth,
		int *dstheight
	);

	/* 

	Zooming functions

	*/

	bitmap_t * transform_zoom(bitmap_t * src, double zoomx, double zoomy, int smooth);

	void get_zoom_size(int width, int height, double zoomx, double zoomy, int *dstwidth, int *dstheight);

	/* 

	Shrinking functions

	*/     

	bitmap_t * transform_shrink(bitmap_t * src, int factorx, int factory);

	/* 

	Specialized rotation functions

	*/

	bitmap_t * transform_rotate_90(bitmap_t * src, int numClockwiseTurns);

	/* Ends C function definitions when using C++ */
#ifdef __cplusplus
}
#endif

#endif				/* LUTRO_ROTOZOOM_H */
