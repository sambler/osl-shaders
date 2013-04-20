/*
 * patterns_osl.h
 * from https://github.com/sambler/osl-shaders
 */

/************************************************************************
 * patterns_osl.h - OSL version of patterns.h (see below)
 * (inc filterwidth.h) for use with Blender/Cycles
 * uses BUILTIN functions where possible
 * removed extern and uniform keywords,
 * change "comp" to index e.g. foo[0],
 * added check for def CCL_STDOSL_H
 * added fadeout function from filterwidth.h
 * filterwidth is now a BUILTIN function
 *
 ************************************************************************/

/************************************************************************
 * patterns.h - Some handy functions for various patterns.  Wherever
 *              possible, antialiased versions will also be given.
 *
 * Author: Larry Gritz (lg AT larrygritz DOT com)
 *
 * Reference:
 *   _Advanced RenderMan: Creating CGI for Motion Picture_,
 *   by Anthony A. Apodaca and Larry Gritz, Morgan Kaufmann, 1999.
 *
 * $Revision: 1.2 $    $Date: 2003/12/24 06:18:06 $
 *
 ************************************************************************/

#ifndef CCL_STDOSL_H
#include "stdosl.h"
#endif

#ifndef PATTERNS_H
#define PATTERNS_H 1

/* fadeout taken from filterwidth.h */
float fadeout(float g, float g_avg, float featuresize, float fwidth) {
        return mix (g, g_avg, smoothstep(.2,.6,fwidth/featuresize))        ;
        }


/* Handy square routine */
float sqr (float x)
{
    return x*x;
}



/* Antialiased abs().
 * Compute the box filter of abs(t) from x-dx/2 to x+dx/2.
 * Hinges on the realization that the indefinite integral of abs(x) is
 * sign(x) * 1/2 x*x;
 */
float filteredabs (float x, float dx)
{
    float integral (float t) {
    return sign(t) * 0.5 * t*t;
    }

    float x0 = x - 0.5*dx;
    float x1 = x0 + dx;
    return (integral(x1) - integral(x0)) / dx;
}




/* Antialiased smoothstep(edge0,edge1,x).
 * Compute the box filter of smoothstep(edge0,edge1,t) from x-dx/2 to x+dx/2.
 * Strategy: divide domain into 3 regions: t < edge0, edge0 <= t <= edge1,
 * and t > edge1.  Region 1 has integral 0.  Region 2 is computed by
 * analytically integrating smoothstep, which is -2t^3+3t^2.  Region 3
 * is trivially 1.
 */
float filteredsmoothstep (float edge0, float edge1, float x, float dx)
{
    float integral (float t) {
        return -0.5*t*t * (t*t - 2*t);
    }

    /* Compute x0, x1 bounding region of integration, and normalize so that
     * edge0==0, edge1==1
     */
    float edgediff = edge1 - edge0;
    float x0 = (x-edge0)/edgediff;
    float fw = dx / edgediff;
    x0 -= 0.5*fw;
    float x1 = x0 + fw;

    /* Region 1 always contributes nothing */
    float integralv = 0;
    /* Region 2 - compute integral in region between 0 and 1 */
    if (x0 < 1 && x1 > 0)
    integralv += integral(min(x1,1)) - integral(max(x0,0));
    /* Region 3 - is 1.0 */
    if (x1 > 1)
    integralv += x1-max(1,x0);
    return integralv / fw;
}



/* A 1-D pulse pattern:  return 1 if edge0 <= x <= edge1, otherwise 0 */
float pulse (float edge0, float edge1, float x)
{
    return step(edge0,x) - step(edge1,x);
}



float filteredpulse (float edge0, float edge1, float x, float dx)
{
    float x0 = x - dx/2;
    float x1 = x0 + dx;
    return max (0, (min(x1,edge1)-max(x0,edge0)) / dx);
}



/* A pulse train: a signal that repeats with a given period, and is
 * 0 when 0 <= mod(x,period) < edge, and 1 when mod(x,period) > edge.
 */
float pulsetrain (float edge, float period, float x)
{
    return pulse (edge, period, mod(x,period));
}


/* Filtered pulse train: it's not as simple as just returning the mod
 * of filteredpulse -- you have to take into account that the filter may
 * cover multiple pulses in the train.
 * Strategy: consider the function that is the integral of the pulse
 * train from 0 to x. Just subtract!
 */
float filteredpulsetrain (float edge, float period, float x, float dx)
{
    /* First, normalize so period == 1 and our domain of interest is > 0 */
    float w = dx/period;
    float x0 = x/period - w/2;
    float x1 = x0+w;
    float nedge = edge / period;   /* normalized edge value */

    /* Definite integral of normalized pulsetrain from 0 to t */
    float integral (float t) {
        float nedge;
        return ((1-nedge)*floor(t) + max(0,t-floor(t)-nedge));
    }

    /* Now we want to integrate the normalized pulsetrain over [x0,x1] */
    return (integral(x1) - integral(x0)) / w;
}



float smoothpulse (float edge0, float edge1, float edge2, float edge3, float x)
{
    return smoothstep(edge0,edge1,x) - smoothstep(edge2,edge3,x);
}


float filteredsmoothpulse (float edge0, float edge1, float edge2, float edge3, float x, float dx)
{
    return filteredsmoothstep(edge0,edge1,x,dx) - filteredsmoothstep(edge2,edge3,x,dx);
}



/* A pulse train of smoothsteps: a signal that repeats with a given
 * period, and is 0 when 0 <= mod(x/period,1) < edge, and 1 when
 * mod(x/period,1) > edge.
 */
float smoothpulsetrain (float edge0, float edge1, float edge2, float edge3, float period, float x)
{
    return smoothpulse (edge0, edge1, edge2, edge3, mod(x,period));
}



/* varyEach takes a computed color, then tweaks each indexed item
 * separately to add some variation.  Hue, saturation, and lightness
 * are all independently controlled.  Hue adds, but saturation and
 * lightness multiply.
 */
color varyEach(color Cin, float index, float varyhue, float varysat, float varylum)
{
    /* Convert to "hsl" space, it's more convenient */
    color Chsl = transformc("hsl", Cin);
    float h = Chsl[0], s = Chsl[1], l = Chsl[2];
    /* Modify Chsl by adding Cvary scaled by our separate h,s,l controls */
    h += varyhue * (cellnoise(index+3)-0.5);
    s *= 1 - varysat * (cellnoise(index-14)-0.5);
    l *= 1 - varylum * (cellnoise(index+37)-0.5);
    Chsl = color (mod(h,1), clamp(s,0,1), clamp(l,0,1));
    /* Clamp hsl and transform back to rgb space */
    return transformc("hsl","rgb", Chsl) ;
}



/* Given 2-D texture coordinates ss,tt and their filter widths ds, dt,
 * and the width and height of the grooves between tiles (assuming that
 * tile spacing is 1.0), figure out which (integer indexed) tile we are
 * on and what coordinates (on [0,1]) within our individual tile we are
 * shading.
 */
float
tilepattern (float ss, float tt, float ds, float dt,
         float groovewidth, float grooveheight,
         output float swhichtile, output float twhichtile,
         output float stile, output float ttile)
{
    swhichtile = floor (ss);
    twhichtile = floor (tt);
    stile = ss - swhichtile;
    ttile = tt - twhichtile;

    return filteredpulsetrain (groovewidth, 1, ss+groovewidth/2, ds)
             * filteredpulsetrain (grooveheight, 1, tt+grooveheight/2, dt);
}



/* basic brick tiling pattern --
 *   inputs:
 *      x, y                    positions on a 2-D surface
 *      tilewidth, tileheight   dimensions of each tile
 *      rowstagger              how much does each row stagger relative to
 *                                   the previous row
 *      rowstaggervary          how much should rowstagger randomly vary
 *      jaggedfreq, jaggedamp   adds noise to the edge between the tiles
 *   outputs:
 *      row, column             index which tile the sample is in
 *      xtile, ytile            position within this tile (0-1)
 */
void basicbrick (float x, float y,
        float tilewidth, float tileheight,
        float rowstagger, float rowstaggervary,
        float jaggedfreq, float jaggedamp,
        output float column, output float row,
        output float xtile, output float ytile
    )
{
    point PP;
    float scoord = x, tcoord = y;

    if (jaggedamp != 0.0) {
    /* Make the shapes of the bricks vary just a bit */
    PP = point(noise(x*jaggedfreq/tilewidth, y*jaggedfreq/tileheight));
    scoord += jaggedamp * PP[0];
    tcoord += jaggedamp * PP[1];
    }

    xtile = scoord / tilewidth;
    ytile = tcoord / tileheight;
    row = floor (ytile);   /* which brick row? */

    /* Shift the columns randomly by row */
    xtile += mod (rowstagger * row, 1);
    xtile += rowstaggervary * (noise (row+0.5) - 0.5);

    column = floor (xtile);
    xtile -= column;
    ytile -= row;
}



#endif /* defined(PATTERNS_H) */

