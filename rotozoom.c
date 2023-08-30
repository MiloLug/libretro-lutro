#include <stdlib.h>
#include <string.h>

#include "rotozoom.h"

/* ---- Internally used structures */

// bytes per pixel
#define BPP 4


typedef struct tColorRGBA {
	uint8_t r;
	uint8_t g;
	uint8_t b;
	uint8_t a;
} tColorRGBA;


#define MAX(a,b) (((a) > (b)) ? (a) : (b))

// Number of guard rows added to destination surfaces to prevent OOB access
#define GUARD_ROWS (2)

// Lower limit of absolute zoom factor or rotation degrees.
#define VALUE_LIMIT	0.00001


/*! 
\brief Internal integer-factor averaging Shrinker.

\param src The surface to shrink (input).
\param dst The shrunken surface (output).
\param factorx The horizontal shrinking ratio.
\param factory The vertical shrinking ratio.

\return 0 for success or -1 for error.
*/
int _transform_shrinkRGBA(bitmap_t * src, bitmap_t * dst, int factorx, int factory)
{
	int x, y, dx, dy, dgap, ra, ga, ba, aa;
	int n_average;
	tColorRGBA *sp, *osp, *oosp;
	tColorRGBA *dp;

	/*
	* Averaging integer shrink
	*/

	/* Precalculate division factor */
	n_average = factorx * factory;

	/*
	* Scan destination
	*/
	sp = (tColorRGBA *) src->data;
	
	dp = (tColorRGBA *) dst->data;
	dgap = dst->pitch - dst->width * 4;

	for (y = 0; y < dst->height; y++) {

		osp=sp;
		for (x = 0; x < dst->width; x++) {

			/* Trace out source box and accumulate */
			oosp=sp;
			ra=ga=ba=aa=0;
			for (dy=0; dy < factory; dy++) {
				for (dx=0; dx < factorx; dx++) {
					ra += sp->r;
					ga += sp->g;
					ba += sp->b;
					aa += sp->a;

					sp++;
				} 
				/* src dx loop */
				sp = (tColorRGBA *)((uint8_t*)sp + (src->pitch - 4*factorx)); // next y
			}
			/* src dy loop */

			/* next box-x */
			sp = (tColorRGBA *)((uint8_t*)oosp + 4*factorx);

			/* Store result in destination */
			dp->r = ra / n_average;
			dp->g = ga / n_average;
			dp->b = ba / n_average;
			dp->a = aa / n_average;

			/*
			* Advance destination pointer 
			*/
			dp++;
		} 
		/* dst x loop */

		/* next box-y */
		sp = (tColorRGBA *)((uint8_t*)osp + src->pitch*factory);

		/*
		* Advance destination pointers 
		*/
		dp = (tColorRGBA *) ((uint8_t *) dp + dgap);
	} 
	/* dst y loop */

	return (0);
}


/*! 
\brief Internal zoomer.

\param src The surface to zoom (input).
\param dst The zoomed surface (output).
\param flipx Flag indicating if the image should be horizontally flipped.
\param flipy Flag indicating if the image should be vertically flipped.
\param smooth Antialiasing flag; set to SMOOTHING_ON to enable.

\return 0 for success or -1 for error.
*/
int _transform_zoomRGBA(const bitmap_t * src, bitmap_t * dst, int flipx, int flipy, bool smooth)
{
	int x, y, sx, sy, ssx, ssy, *sax, *say, *csax, *csay, *salast, csx, csy, ex, ey, cx, cy, sstep, sstepx, sstepy;
	tColorRGBA *c00, *c01, *c10, *c11;
	tColorRGBA *sp, *csp, *dp;
	int spixelgap, spixelw, spixelh, dgap, t1, t2;

	/*
	* Allocate memory for row/column increments 
	*/
	if ((sax = (int *) malloc((dst->width + 1) * sizeof(uint32_t))) == NULL) {
		return (-1);
	}
	if ((say = (int *) malloc((dst->height + 1) * sizeof(uint32_t))) == NULL) {
		free(sax);
		return (-1);
	}

	/*
	* Precalculate row increments 
	*/
	spixelw = (src->width - 1);
	spixelh = (src->height - 1);
	if (smooth) {
		sx = (int) (65536.0 * (float) spixelw / (float) (dst->width - 1));
		sy = (int) (65536.0 * (float) spixelh / (float) (dst->height - 1));
	} else {
		sx = (int) (65536.0 * (float) (src->width) / (float) (dst->width));
		sy = (int) (65536.0 * (float) (src->height) / (float) (dst->height));
	}

	/* Maximum scaled source size */
	ssx = (src->width << 16) - 1;
	ssy = (src->height << 16) - 1;

	/* Precalculate horizontal row increments */
	csx = 0;
	csax = sax;
	for (x = 0; x <= dst->width; x++) {
		*csax = csx;
		csax++;
		csx += sx;

		/* Guard from overflows */
		if (csx > ssx) { 
			csx = ssx; 
		}
	}

	/* Precalculate vertical row increments */
	csy = 0;
	csay = say;
	for (y = 0; y <= dst->height; y++) {
		*csay = csy;
		csay++;
		csy += sy;

		/* Guard from overflows */
		if (csy > ssy) {
			csy = ssy;
		}
	}

	sp = (tColorRGBA *) src->data;
	dp = (tColorRGBA *) dst->data;
	dgap = dst->pitch - dst->width * 4;
	spixelgap = src->pitch/4;

	if (flipx) sp += spixelw;
	if (flipy) sp += (spixelgap * spixelh);

	/*
	* Switch between interpolating and non-interpolating code 
	*/
	if (smooth) {

		/*
		* Interpolating Zoom 
		*/
		csay = say;
		for (y = 0; y < dst->height; y++) {
			csp = sp;
			csax = sax;
			for (x = 0; x < dst->width; x++) {
				/*
				* Setup color source pointers 
				*/
				ex = (*csax & 0xffff);
				ey = (*csay & 0xffff);
				cx = (*csax >> 16);
				cy = (*csay >> 16);
				sstepx = cx < spixelw;
				sstepy = cy < spixelh;
				c00 = sp;
				c01 = sp;
				c10 = sp;
				if (sstepy) {
					if (flipy) {
						c10 -= spixelgap;
					} else {
						c10 += spixelgap;
					}
				}
				c11 = c10;
				if (sstepx) {
					if (flipx) {
						c01--;
						c11--;
					} else {
						c01++;
						c11++;
					}
				}

				/*
				* Draw and interpolate colors 
				*/
				t1 = ((((c01->r - c00->r) * ex) >> 16) + c00->r) & 0xff;
				t2 = ((((c11->r - c10->r) * ex) >> 16) + c10->r) & 0xff;
				dp->r = (((t2 - t1) * ey) >> 16) + t1;
				t1 = ((((c01->g - c00->g) * ex) >> 16) + c00->g) & 0xff;
				t2 = ((((c11->g - c10->g) * ex) >> 16) + c10->g) & 0xff;
				dp->g = (((t2 - t1) * ey) >> 16) + t1;
				t1 = ((((c01->b - c00->b) * ex) >> 16) + c00->b) & 0xff;
				t2 = ((((c11->b - c10->b) * ex) >> 16) + c10->b) & 0xff;
				dp->b = (((t2 - t1) * ey) >> 16) + t1;
				t1 = ((((c01->a - c00->a) * ex) >> 16) + c00->a) & 0xff;
				t2 = ((((c11->a - c10->a) * ex) >> 16) + c10->a) & 0xff;
				dp->a = (((t2 - t1) * ey) >> 16) + t1;				
				/*
				* Advance source pointer x
				*/
				salast = csax;
				csax++;				
				sstep = (*csax >> 16) - (*salast >> 16);
				if (flipx) {
					sp -= sstep;
				} else {
					sp += sstep;
				}

				/*
				* Advance destination pointer x
				*/
				dp++;
			}
			/*
			* Advance source pointer y
			*/
			salast = csay;
			csay++;
			sstep = (*csay >> 16) - (*salast >> 16);
			sstep *= spixelgap;
			if (flipy) { 
				sp = csp - sstep;
			} else {
				sp = csp + sstep;
			}

			/*
			* Advance destination pointer y
			*/
			dp = (tColorRGBA *) ((uint8_t *) dp + dgap);
		}
	} else {
		/*
		* Non-Interpolating Zoom 
		*/		
		csay = say;
		for (y = 0; y < dst->height; y++) {
			csp = sp;
			csax = sax;
			for (x = 0; x < dst->width; x++) {
				/*
				* Draw 
				*/
				*dp = *sp;

				/*
				* Advance source pointer x
				*/
				salast = csax;
				csax++;				
				sstep = (*csax >> 16) - (*salast >> 16);
				if (flipx) sstep = -sstep;
				sp += sstep;

				/*
				* Advance destination pointer x
				*/
				dp++;
			}
			/*
			* Advance source pointer y
			*/
			salast = csay;
			csay++;
			sstep = (*csay >> 16) - (*salast >> 16);
			sstep *= spixelgap;
			if (flipy) sstep = -sstep;			
			sp = csp + sstep;

			/*
			* Advance destination pointer y
			*/
			dp = (tColorRGBA *) ((uint8_t *) dp + dgap);
		}
	}

	/*
	* Remove temp arrays 
	*/
	free(sax);
	free(say);

	return (0);
}

/*! 
\brief Internal rotozoomer with optional anti-aliasing.

\param src Source surface.
\param dst Destination surface.
\param cx Horizontal center coordinate.
\param cy Vertical center coordinate.
\param isin Integer version of sine of angle.
\param icos Integer version of cosine of angle.
\param flipx Flag indicating horizontal mirroring should be applied.
\param flipy Flag indicating vertical mirroring should be applied.
\param smooth Flag indicating anti-aliasing should be used.
*/
void _transformSurfaceRGBA(const bitmap_t * src, bitmap_t * dst, int32_t cx, int32_t cy, int32_t isin, int32_t icos, bool flipx, bool flipy, bool smooth)
{
	int32_t x, y, t1, t2, dx, dy, xd, yd, sdx, sdy, ax, ay, ex, ey, sw, sh;
	tColorRGBA c00, c01, c10, c11, cswap;
	tColorRGBA *pc, *sp;
	int32_t gap;

	/*
	* Variable setup 
	*/
	xd = ((src->width - dst->width) << 15);
	yd = ((src->height - dst->height) << 15);
	ax = (cx << 16) - (icos * cx);
	ay = (cy << 16) - (isin * cx);
	sw = src->width - 1;
	sh = src->height - 1;
	pc = (tColorRGBA*) dst->data;
	gap = dst->pitch - dst->width * 4;

	/*
	* Switch between interpolating and non-interpolating code 
	*/
	if (smooth) {
		for (y = 0; y < dst->height; y++) {
			dy = cy - y;
			sdx = (ax + (isin * dy)) + xd;
			sdy = (ay - (icos * dy)) + yd;
			for (x = 0; x < dst->width; x++) {
				dx = (sdx >> 16);
				dy = (sdy >> 16);
				if (flipx) dx = sw - dx;
				if (flipy) dy = sh - dy;
				if ((dx > -1) && (dy > -1) && (dx < (src->width-1)) && (dy < (src->height-1))) {
					sp = (tColorRGBA *)src->data;
					sp += ((src->pitch/4) * dy);
					sp += dx;
					c00 = *sp;
					sp += 1;
					c01 = *sp;
					sp += (src->pitch/4);
					c11 = *sp;
					sp -= 1;
					c10 = *sp;
					if (flipx) {
						cswap = c00; c00=c01; c01=cswap;
						cswap = c10; c10=c11; c11=cswap;
					}
					if (flipy) {
						cswap = c00; c00=c10; c10=cswap;
						cswap = c01; c01=c11; c11=cswap;
					}
					/*
					* Interpolate colors 
					*/
					ex = (sdx & 0xffff);
					ey = (sdy & 0xffff);
					t1 = ((((c01.r - c00.r) * ex) >> 16) + c00.r) & 0xff;
					t2 = ((((c11.r - c10.r) * ex) >> 16) + c10.r) & 0xff;
					pc->r = (((t2 - t1) * ey) >> 16) + t1;
					t1 = ((((c01.g - c00.g) * ex) >> 16) + c00.g) & 0xff;
					t2 = ((((c11.g - c10.g) * ex) >> 16) + c10.g) & 0xff;
					pc->g = (((t2 - t1) * ey) >> 16) + t1;
					t1 = ((((c01.b - c00.b) * ex) >> 16) + c00.b) & 0xff;
					t2 = ((((c11.b - c10.b) * ex) >> 16) + c10.b) & 0xff;
					pc->b = (((t2 - t1) * ey) >> 16) + t1;
					t1 = ((((c01.a - c00.a) * ex) >> 16) + c00.a) & 0xff;
					t2 = ((((c11.a - c10.a) * ex) >> 16) + c10.a) & 0xff;
					pc->a = (((t2 - t1) * ey) >> 16) + t1;
				}
				sdx += icos;
				sdy += isin;
				pc++;
			}
			pc = (tColorRGBA *) ((uint8_t *) pc + gap);
		}
	} else {
		for (y = 0; y < dst->height; y++) {
			dy = cy - y;
			sdx = (ax + (isin * dy)) + xd;
			sdy = (ay - (icos * dy)) + yd;
			for (x = 0; x < dst->width; x++) {
				dx = (short) (sdx >> 16);
				dy = (short) (sdy >> 16);
				if (flipx) dx = (src->width-1)-dx;
				if (flipy) dy = (src->height-1)-dy;
				if ((dx >= 0) && (dy >= 0) && (dx < src->width) && (dy < src->height)) {
					sp = (tColorRGBA *) ((uint8_t *) src->data + src->pitch * dy);
					sp += dx;
					*pc = *sp;
				}
				sdx += icos;
				sdy += isin;
				pc++;
			}
			pc = (tColorRGBA *) ((uint8_t *) pc + gap);
		}
	}
}

/*!
\brief Rotates in increments of 90 degrees.

\param src Source surface to rotate.
\param numClockwiseTurns Number of clockwise 90 degree turns to apply to the source.

\returns The new, rotated surface; or NULL for surfaces with incorrect input format.
*/
bitmap_t * transform_rotate_90(bitmap_t * src, int numClockwiseTurns) 
{
	int row, col, newWidth, newHeight;
	int bpr;
	bitmap_t * dst;
	uint8_t* srcBuf;
	uint8_t* dstBuf;
	int normalizedClockwiseTurns;

	if (!src)
	    return NULL;

	/* normalize numClockwiseTurns */
	normalizedClockwiseTurns = (numClockwiseTurns % 4);
	if (normalizedClockwiseTurns < 0) {
		normalizedClockwiseTurns += 4;
	}

	/* If turns are even, our new width/height will be the same as the source surface */
	if (normalizedClockwiseTurns % 2) {
		newWidth = src->height;
		newHeight = src->width;
	} else {
		newWidth = src->width;
		newHeight = src->height;
	}

	dst = pntr_new_bmp(newWidth, newHeight);

	switch(normalizedClockwiseTurns) {
		case 0: /* Make a copy of the surface */
			{
				if (src->pitch == dst->pitch) {
					/* If the pitch is the same for both surfaces, the memory can be copied all at once. */
					memcpy(dst->data, src->data, (src->height * src->pitch));
				}
				else
				{
					/* If the pitch differs, copy each row separately */
					srcBuf = (uint8_t*)(src->data);
					dstBuf = (uint8_t*)(dst->data);
					bpr = src->width * BPP;
					for (row = 0; row < src->height; row++) {
						memcpy(dstBuf, srcBuf, bpr);
						srcBuf += src->pitch;
						dstBuf += dst->pitch;
					}
				}
			}
			break;

			/* rotate clockwise */
		case 1: /* rotated 90 degrees clockwise */
			{
				for (row = 0; row < src->height; ++row) {
					srcBuf = (uint8_t*)(src->data) + (row * src->pitch);
					dstBuf = (uint8_t*)(dst->data) + (dst->width - row - 1) * BPP;
					for (col = 0; col < src->width; ++col) {
						memcpy (dstBuf, srcBuf, BPP);
						srcBuf += BPP;
						dstBuf += dst->pitch;
					} 
				} 
			}
			break;

		case 2: /* rotated 180 degrees clockwise */
			{
				for (row = 0; row < src->height; ++row) {
					srcBuf = (uint8_t*)(src->data) + (row * src->pitch);
					dstBuf = (uint8_t*)(dst->data) + ((dst->height - row - 1) * dst->pitch) + (dst->width - 1) * BPP;
					for (col = 0; col < src->width; ++col) {
						memcpy (dstBuf, srcBuf, BPP);
						srcBuf += BPP;
						dstBuf -= BPP;
					} 
				} 
			}
			break;

		case 3: /* rotated 270 degrees clockwise */
			{
				for (row = 0; row < src->height; ++row) {
					srcBuf = (uint8_t*)(src->data) + (row * src->pitch);
					dstBuf = (uint8_t*)(dst->data) + (row * BPP) + ((dst->height - 1) * dst->pitch);
					for (col = 0; col < src->width; ++col) {
						memcpy (dstBuf, srcBuf, BPP);
						srcBuf += BPP;
						dstBuf -= dst->pitch;
					} 
				} 
			}
			break;
	}

	return dst;
}


/*!
\brief Calculates the size of the target as if rotozoom was already applied. 

\param width The source surface width.
\param height The source surface height.
\param angle The angle to rotate in radianss.
\param zoomx The horizontal scaling factor.
\param zoomy The vertical scaling factor.
\param dstwidth The calculated width of the destination surface.
\param dstheight The calculated height of the destination surface.
\param canglezoom The sine of the angle adjusted by the zoom factor.
\param sanglezoom The cosine of the angle adjusted by the zoom factor.

*/
void _get_rotozoom_size_trig(int width, int height, double angle, double zoomx, double zoomy, 
	int *dstwidth, int *dstheight, 
	double *canglezoom, double *sanglezoom)
{
	double x, y, cx, cy, sx, sy;
	int dstwidthhalf, dstheighthalf;

	/*
	* Determine destination width and height by rotating a centered source box 
	*/
	*sanglezoom = sin(angle);
	*canglezoom = cos(angle);
	*sanglezoom *= zoomx;
	*canglezoom *= zoomy;
	x = (double)(width / 2);
	y = (double)(height / 2);
	cx = *canglezoom * x;
	cy = *canglezoom * y;
	sx = *sanglezoom * x;
	sy = *sanglezoom * y;

	dstwidthhalf = MAX((int)
		ceil(MAX(MAX(MAX(fabs(cx + sy), fabs(cx - sy)), fabs(-cx + sy)), fabs(-cx - sy))), 1);
	dstheighthalf = MAX((int)
		ceil(MAX(MAX(MAX(fabs(sx + cy), fabs(sx - cy)), fabs(-sx + cy)), fabs(-sx - cy))), 1);
	*dstwidth = 2 * dstwidthhalf;
	*dstheight = 2 * dstheighthalf;
}

void get_rotozoom_size_xy(int width, int height, double angle, double zoomx, double zoomy, int *dstwidth, int *dstheight)
{
	double dummy_sanglezoom, dummy_canglezoom;

	_get_rotozoom_size_trig(width, height, angle, zoomx, zoomy, dstwidth, dstheight, &dummy_sanglezoom, &dummy_canglezoom);
}

void get_rotozoom_size(int width, int height, double angle, double zoom, int *dstwidth, int *dstheight)
{
	double dummy_sanglezoom, dummy_canglezoom;

	_get_rotozoom_size_trig(width, height, angle, zoom, zoom, dstwidth, dstheight, &dummy_sanglezoom, &dummy_canglezoom);
}

/*!
\brief Rotates and zooms. 

\param src The surface to rotozoom.
\param angle The angle to rotate in radians.
\param zoom The scaling factor.
\param smooth Antialiasing flag; set to SMOOTHING_ON to enable.

\return The new rotozoomed surface.
*/
bitmap_t *transform_rotozoom(const bitmap_t * src, double angle, double zoom, bool smooth)
{
	return transform_rotozoom_xy(src, angle, zoom, zoom, smooth);
}

/*!
\brief Rotates and zooms with different horizontal and vertival scaling factors and optional anti-aliasing. 

\param src The surface to rotozoom.
\param angle The angle to rotate in radians.
\param zoomx The horizontal scaling factor.
\param zoomy The vertical scaling factor.
\param smooth Antialiasing flag; set to SMOOTHING_ON to enable.

\return The new rotozoomed surface.
*/
bitmap_t *transform_rotozoom_xy(const bitmap_t * src, double angle, double zoomx, double zoomy, bool smooth)
{
	const bitmap_t *rz_src = src;
	bitmap_t *rz_dst;
	double zoominv;
	double sanglezoom, canglezoom, sanglezoominv, canglezoominv;
	int32_t dstwidthhalf, dstwidth, dstheighthalf, dstheight;
	bool flipx, flipy;

	/*
	* Sanity check 
	*/
	if (src == NULL)
		return NULL;

	/*
	* Sanity check zoom factor 
	*/
	flipx = zoomx < 0.0;
	if (flipx) zoomx = -zoomx;
	flipy = zoomy < 0.0;
	if (flipy) zoomy = -zoomy;
	if (zoomx < VALUE_LIMIT) zoomx = VALUE_LIMIT;
	if (zoomy < VALUE_LIMIT) zoomy = VALUE_LIMIT;
	zoominv = 65536.0 / (zoomx * zoomx);

	/*
	* Check if we have a rotozoom or just a zoom 
	*/
	if (fabs(angle) > VALUE_LIMIT) {
		/*
		* Angle!=0: full rotozoom 
		*/

		_get_rotozoom_size_trig(rz_src->width, rz_src->height, angle, zoomx, zoomy, &dstwidth, &dstheight, &canglezoom, &sanglezoom);

		/*
		* Calculate target factors from sin/cos and zoom 
		*/
		sanglezoominv = sanglezoom;
		canglezoominv = canglezoom;
		sanglezoominv *= zoominv;
		canglezoominv *= zoominv;

		/* Calculate half size */
		dstwidthhalf = dstwidth / 2;
		dstheighthalf = dstheight / 2;

		rz_dst = pntr_new_bmp(dstwidth, dstheight + GUARD_ROWS);
		rz_dst->height = dstheight;  // Adjust for guard rows

		_transformSurfaceRGBA(
			rz_src,
			rz_dst,
			dstwidthhalf,
			dstheighthalf,
			(int) (sanglezoominv),
			(int) (canglezoominv), 
			flipx,
			flipy,
			smooth
		);
	} else {
		/*
		* Angle=0: Just a zoom 
		*/
		
		get_zoom_size(rz_src->width, rz_src->height, zoomx, zoomy, &dstwidth, &dstheight);

		rz_dst = pntr_new_bmp(dstwidth, dstheight + GUARD_ROWS);
		rz_dst->height = dstheight;  // Adjust for guard rows

		_transform_zoomRGBA(rz_src, rz_dst, flipx, flipy, smooth);
	}

	return rz_dst;
}

/*!
\brief Calculates the size of the target as if zoom was already applied.

\param width The width of the source surface to zoom.
\param height The height of the source surface to zoom.
\param zoomx The horizontal zoom factor.
\param zoomy The vertical zoom factor.
\param dstwidth Pointer to an integer to store the calculated width of the zoomed target surface.
\param dstheight Pointer to an integer to store the calculated height of the zoomed target surface.
*/
void get_zoom_size(int width, int height, double zoomx, double zoomy, int *dstwidth, int *dstheight)
{
	/*
	* Make zoom factors positive 
	*/
	int flipx, flipy;
	flipx = (zoomx<0.0);
	if (flipx) zoomx = -zoomx;
	flipy = (zoomy<0.0);
	if (flipy) zoomy = -zoomy;

	/*
	* Sanity check zoom factors 
	*/
	if (zoomx < VALUE_LIMIT) {
		zoomx = VALUE_LIMIT;
	}
	if (zoomy < VALUE_LIMIT) {
		zoomy = VALUE_LIMIT;
	}

	/*
	* Calculate target size 
	*/
	*dstwidth = (int) floor(((double) width * zoomx) + 0.5);
	*dstheight = (int) floor(((double) height * zoomy) + 0.5);
	if (*dstwidth < 1) {
		*dstwidth = 1;
	}
	if (*dstheight < 1) {
		*dstheight = 1;
	}
}

/*! 
\brief Zoom by independent horizontal and vertical factors with optional smoothing.

\param src The surface to zoom.
\param zoomx The horizontal zoom factor.
\param zoomy The vertical zoom factor.
\param smooth Antialiasing flag; set to SMOOTHING_ON to enable.

\return The new, zoomed surface.
*/
bitmap_t *transform_zoom(bitmap_t * src, double zoomx, double zoomy, bool smooth)
{
	bitmap_t *rz_src = src;
	bitmap_t *rz_dst;
	int dstwidth, dstheight;
	int flipx, flipy;

	/*
	* Sanity check 
	*/
	if (src == NULL)
		return (NULL);

	flipx = zoomx < 0.0;
	if (flipx) zoomx = -zoomx;
	flipy = zoomy < 0.0;
	if (flipy) zoomy = -zoomy;

	get_zoom_size(rz_src->width, rz_src->height, zoomx, zoomy, &dstwidth, &dstheight);

	rz_dst = pntr_new_bmp(dstwidth, dstheight + GUARD_ROWS);
	rz_dst->height = dstheight;  // Adjust for guard rows

	_transform_zoomRGBA(rz_src, rz_dst, flipx, flipy, smooth);

	return rz_dst;
}

/*! 
\brief Shrink by an integer ratio using averaging.

'factorx' and 'factory' are the shrinking ratios (i.e. 2=1/2 the size,
3=1/3 the size, etc.)

\param src The surface to shrink.
\param factorx The horizontal shrinking ratio.
\param factory The vertical shrinking ratio.

\return The new, shrunken surface.
*/
/*@null@*/ 
bitmap_t *transform_shrink(bitmap_t *src, int factorx, int factory)
{
	int result;
	bitmap_t *rz_src = src;
	bitmap_t *rz_dst = NULL;
	int dstwidth, dstheight;

	if (src == NULL)
		return NULL;

	/* Get size for target */
	dstwidth=rz_src->width / factorx;
	while (dstwidth*factorx > rz_src->width) dstwidth--;
	dstheight=rz_src->height / factory;
	while (dstheight*factory > rz_src->height) dstheight--;

	rz_dst = pntr_new_bmp(dstwidth, dstheight + GUARD_ROWS);
	rz_dst->height = dstheight;  // Adjust for guard rows

	result = _transform_shrinkRGBA(rz_src, rz_dst, factorx, factory);		
	if ((result!=0) || (rz_dst==NULL))
		rz_dst = NULL;

	return rz_dst;
}
