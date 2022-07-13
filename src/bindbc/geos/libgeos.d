module bindbc.geos.libgeos;

/**
 * Custom added part to support multiple versions of the library
 */

enum GEOSSupport
{
    geos_3_8, // lowest supported version
    geos_3_10
}

version (GEOS_3_8) enum geosSupport = GEOSSupport.geos_3_8;
else version (GEOS_3_10) enum geosSupport = GEOSSupport.geos_3_10;
else enum geosSupport = GEOSSupport.geos_3_8;

nothrow @nogc:

/************************************************************************
 *
 * C-Wrapper for GEOS library
 *
 * Copyright (C) 2010 2011 Sandro Santilli <strk@kbt.io>
 * Copyright (C) 2005 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 ***********************************************************************/

/**
* \file geos_c.h
* \brief C API for the GEOS geometry algorithms library.
*
* The C API is the preferred API to use when integration GEOS into
* you program/language/etc. While the C++ API is available, the ABI
* will not be stable between versions, and the API may also change.
* The GEOS team makes an effort to keep the C API stable, and to
* deprecate function signatures only over a long time period to allow
* transition time.
*
* Important programming notes:
*
* - Remember to call initGEOS() before any use of this library's
*   functions, and call finishGEOS() when done.
* - Currently you have to explicitly GEOSGeom_destroy() all
*   GEOSGeom objects to avoid memory leaks, and GEOSFree()
*   all returned char * (unless const).
* - Functions ending with _r are thread safe (reentrant);
*   see details in http://trac.osgeo.org/geos/wiki/RFC3.
*   To avoid accidental use of non-reentrant functions,
*   define GEOS_USE_ONLY_R_API before including geos_c.h.
*
*/

extern (C):

/* for size_t definition */

/* ====================================================================== */
/* Version */
/* ====================================================================== */

/** \cond */

enum GEOS_VERSION_MAJOR = 3;

enum GEOS_VERSION_MINOR = 10;

enum GEOS_VERSION_PATCH = 2;

enum GEOS_VERSION = "3.10.2";

enum GEOS_JTS_PORT = "1.18.0";

enum GEOS_CAPI_VERSION_MAJOR = 1;
enum GEOS_CAPI_VERSION_MINOR = 16;
enum GEOS_CAPI_VERSION_PATCH = 0;
enum GEOS_CAPI_VERSION = "3.10.2-CAPI-1.16.0";

enum GEOS_CAPI_FIRST_INTERFACE = GEOS_CAPI_VERSION_MAJOR;
enum GEOS_CAPI_LAST_INTERFACE = GEOS_CAPI_VERSION_MAJOR + GEOS_CAPI_VERSION_MINOR;

/** \endcond */

/**
* Type returned by GEOS_init_r(), for use in multi-threaded
* applications.
*
* There should be only one GEOSContextHandle_t per thread.
*/
struct GEOSContextHandle_HS;
alias GEOSContextHandle_t = GEOSContextHandle_HS*;

/**
* Callback function for passing GEOS error messages to parent process.
*
* Set the GEOSMessageHandler for error and notice messages in \ref initGEOS
* for single-threaded programs, or using \ref initGEOS_r for threaded
* programs
*
* \param fmt the message format template
*/
alias GEOSMessageHandler = void function (const(char)* fmt, ...);

/**
* A GEOS message handler function.
*
* \param message the message contents
* \param userdata the user data pointer that was passed to GEOS when
* registering this message handler.
*
* \see GEOSContext_setErrorMessageHandler
* \see GEOSContext_setNoticeMessageHandler
*/
alias GEOSMessageHandler_r = void function (const(char)* message, void* userdata);

/*
* When we're included by geos_c.cpp, these types are #defined to the
* C++ definitions via preprocessor. We don't touch them to allow the
* compiler to cross-check the declarations. However, for all "normal"
* C-API users, we need to define these types as "opaque" struct pointers, as
* those clients don't have access to the original C++ headers, by design.
*/

/**
* Geometry generic type. Geometry can be a point, linestring, polygon,
* multipoint, multilinestring, multipolygon, or geometrycollection.
* Geometry type can be read with \ref GEOSGeomTypeId. Most functions
* in GEOS either have GEOSGeometry* as a parameter or a return type.
* \see GEOSGeom_createPoint
* \see GEOSGeom_createLineString
* \see GEOSGeom_createPolygon
* \see GEOSGeom_createCollection
* \see GEOSGeom_destroy
*/
struct GEOSGeom_t;
alias GEOSGeometry = GEOSGeom_t;

/**
* Prepared geometry type.
* \see GEOSPrepare()
* \see GEOSPreparedGeom_destroy()
*/
struct GEOSPrepGeom_t;
alias GEOSPreparedGeometry = GEOSPrepGeom_t;

/**
* Coordinate sequence.
* \see GEOSCoordSeq_create()
* \see GEOSCoordSeq_destroy()
*/
struct GEOSCoordSeq_t;
alias GEOSCoordSequence = GEOSCoordSeq_t;

/**
* STRTree index.
* \see GEOSSTRtree_create()
* \see GEOSSTRtree_destroy()
*/
struct GEOSSTRtree_t;
alias GEOSSTRtree = GEOSSTRtree_t;

/**
* Parameter object for buffering.
* \see GEOSBufferParams_create()
* \see GEOSBufferParams_destroy()
*/
struct GEOSBufParams_t;
alias GEOSBufferParams = GEOSBufParams_t;

/**
* Parameter object for validity enforcement.
* \see GEOSMakeValidParams_create()
* \see GEOSMakeValidParams_destroy()
*/
struct GEOSMakeValidParams_t;
alias GEOSMakeValidParams = GEOSMakeValidParams_t;

/** \cond */

/*
* These are compatibility definitions for source compatibility
* with GEOS 2.X clients relying on that type.
*/
alias GEOSGeom = GEOSGeom_t*;
alias GEOSCoordSeq = GEOSCoordSeq_t*;

/** \endcond */

/**
* Geometry type number, used by functions returning or
* consuming geometry types.
*
* \see GEOSGeomType
* \see GEOSGeomTypeId
*/
enum GEOSGeomTypes
{
    /** Point */
    GEOS_POINT = 0,
    /** Linestring */
    GEOS_LINESTRING = 1,
    /** Linear ring, used within polygons */
    GEOS_LINEARRING = 2,
    /** Polygon */
    GEOS_POLYGON = 3,
    /** Multipoint, a homogeneous collection of points */
    GEOS_MULTIPOINT = 4,
    /** Multilinestring, a homogeneous collection of linestrings */
    GEOS_MULTILINESTRING = 5,
    /** Multipolygon, a homogeneous collection of polygons */
    GEOS_MULTIPOLYGON = 6,
    /** Geometry collection, a heterogeneous collection of geometry */
    GEOS_GEOMETRYCOLLECTION = 7
}

/**
* Well-known binary byte orders used when
* writing to WKB.
*
* \see GEOSWKBWriter_setByteOrder
*/
enum GEOSWKBByteOrders
{
    /** Big Endian */
    GEOS_WKB_XDR = 0,
    /** Little Endian */
    GEOS_WKB_NDR = 1
}

/**
* Well-known binary flavors to use
* when writing to WKB. ISO flavour is
* more standard. Extended flavour supports
* 3D and SRID embedding. GEOS reads both
* transparently.
*
* \see GEOSWKBWriter_setFlavor
*/
enum GEOSWKBFlavors
{
    /** Extended */
    GEOS_WKB_EXTENDED = 1,
    /** ISO */
    GEOS_WKB_ISO = 2
}

/**
* Callback function for use in spatial index search calls. Pass into
* the query function and handle query results as the index
* returns them.
*
* \see GEOSSTRtree_query
*/
alias GEOSQueryCallback = void function (void* item, void* userdata);

/**
* Callback function for use in spatial index nearest neighbor calculations.
* Allows custom distance to be calculated between items in the
* index. Is passed two items, and sets the calculated distance
* between the items into the distance pointer. Extra data for the
* calculation can be passed via the userdata.
*
* \param item1 first of the pair of items to calculate distance between
* \param item2 second of the pair of items to calculate distance between
* \param distance the distance between the items here
* \param userdata extra data for the calculation
*
* \return zero if distance calculation succeeded, non-zero otherwise
*
* \see GEOSSTRtree_nearest_generic
* \see GEOSSTRtree_iterate
*/
alias GEOSDistanceCallback = int function (
    const(void)* item1,
    const(void)* item2,
    double* distance,
    void* userdata);

/* ========== Interruption ========== */

/**
* Callback function for use in interruption. The callback will be invoked _before_ checking for
* interruption, so can be used to request it.
*
* \see GEOS_interruptRegisterCallback
* \see GEOS_interruptRequest
* \see GEOS_interruptCancel
*/
alias GEOSInterruptCallback = void function ();

/**
* Register a function to be called when processing is interrupted.
* \param cb Callback function to invoke
* \return the previously configured callback
* \see GEOSInterruptCallback
*/
void function (void function () cb) GEOS_interruptRegisterCallback (
    void function () cb);

/**
* Request safe interruption of operations
*/
void GEOS_interruptRequest ();

/**
* Cancel a pending interruption request
*/
void GEOS_interruptCancel ();

/* ========== Initialization and Cleanup ========== */

/**
* Initialize a context for this thread. Pass this context into
* your other calls of `*_r` functions.
* \return a GEOS context for this thread
*/
GEOSContextHandle_t GEOS_init_r ();

/**
* Free the memory associated with a \ref GEOSContextHandle_t
* when you are finished calling GEOS functions.
* \param handle to be freed
*/
void GEOS_finish_r (GEOSContextHandle_t handle);

/**
* Set the notice handler callback function for run-time notice messages.
* \param extHandle the context returned by \ref GEOS_init_r.
* \param nf the handler callback
* \return the previously configured message handler or NULL if no message handler was configured
*/
GEOSMessageHandler GEOSContext_setNoticeHandler_r (
    GEOSContextHandle_t extHandle,
    GEOSMessageHandler nf);

/**  */

/**
* Set the notice handler callback function for run-time error messages.
* \param extHandle the GEOS context from \ref GEOS_init_r
* \param ef the handler callback
* \return the previously configured message handler or NULL if no message handler was configured
*/
GEOSMessageHandler GEOSContext_setErrorHandler_r (
    GEOSContextHandle_t extHandle,
    GEOSMessageHandler ef);

/**
* Sets a notice message handler on the given GEOS context.
* \param extHandle the GEOS context from \ref GEOS_init_r
* \param nf the message handler
* \param userData optional user data pointer that will be passed to the message handler
* \return the previously configured message handler or NULL if no message handler was configured
*/
GEOSMessageHandler_r GEOSContext_setNoticeMessageHandler_r (
    GEOSContextHandle_t extHandle,
    GEOSMessageHandler_r nf,
    void* userData);

/**
* Sets an error message handler on the given GEOS context.
*
* \param extHandle the GEOS context
* \param ef the message handler
* \param userData optional user data pointer that will be passed to the message handler
*
* \return the previously configured message handler or NULL if no message handler was configured
*/
GEOSMessageHandler_r GEOSContext_setErrorMessageHandler_r (
    GEOSContextHandle_t extHandle,
    GEOSMessageHandler_r ef,
    void* userData);

/* ========== Coordinate Sequence functions ========== */

/** \see GEOSCoordSeq_create */
GEOSCoordSequence* GEOSCoordSeq_create_r (
    GEOSContextHandle_t handle,
    uint size,
    uint dims);

/** \see GEOSCoordSeq_copyFromBuffer */
GEOSCoordSequence* GEOSCoordSeq_copyFromBuffer_r (
    GEOSContextHandle_t handle,
    const(double)* buf,
    uint size,
    int hasZ,
    int hasM);

/** \see GEOSCoordSeq_copyFromArrays */
GEOSCoordSequence* GEOSCoordSeq_copyFromArrays_r (
    GEOSContextHandle_t handle,
    const(double)* x,
    const(double)* y,
    const(double)* z,
    const(double)* m,
    uint size);

/** \see GEOSCoordSeq_copyToBuffer */
int GEOSCoordSeq_copyToBuffer_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    double* buf,
    int hasZ,
    int hasM);

/** \see GEOSCoordSeq_copyToArrays */
int GEOSCoordSeq_copyToArrays_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    double* x,
    double* y,
    double* z,
    double* m);

/** \see GEOSCoordSeq_clone */
GEOSCoordSequence* GEOSCoordSeq_clone_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s);

/** \see GEOSCoordSeq_destroy */
void GEOSCoordSeq_destroy_r (GEOSContextHandle_t handle, GEOSCoordSequence* s);

/** \see GEOSCoordSeq_setX */
int GEOSCoordSeq_setX_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s,
    uint idx,
    double val);

/** \see GEOSCoordSeq_setY */
int GEOSCoordSeq_setY_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s,
    uint idx,
    double val);

/** \see GEOSCoordSeq_setZ */
int GEOSCoordSeq_setZ_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s,
    uint idx,
    double val);

/** \see GEOSCoordSeq_setXY */
int GEOSCoordSeq_setXY_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s,
    uint idx,
    double x,
    double y);

/** \see GEOSCoordSeq_setXYZ */
int GEOSCoordSeq_setXYZ_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s,
    uint idx,
    double x,
    double y,
    double z);

/** \see GEOSCoordSeq_setOrdinate */
int GEOSCoordSeq_setOrdinate_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s,
    uint idx,
    uint dim,
    double val);

/** \see GEOSCoordSeq_getX */
int GEOSCoordSeq_getX_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint idx,
    double* val);

/** \see GEOSCoordSeq_getY */
int GEOSCoordSeq_getY_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint idx,
    double* val);

/** \see GEOSCoordSeq_getZ */
int GEOSCoordSeq_getZ_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint idx,
    double* val);

/** \see GEOSCoordSeq_getXY */
int GEOSCoordSeq_getXY_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint idx,
    double* x,
    double* y);

/** \see GEOSCoordSeq_getXYZ */
int GEOSCoordSeq_getXYZ_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint idx,
    double* x,
    double* y,
    double* z);

/** \see GEOSCoordSeq_getOrdinate */
int GEOSCoordSeq_getOrdinate_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint idx,
    uint dim,
    double* val);

/** \see GEOSCoordSeq_getSize */
int GEOSCoordSeq_getSize_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint* size);

/** \see GEOSCoordSeq_getDimensions */
int GEOSCoordSeq_getDimensions_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    uint* dims);

/** \see GEOSCoordSeq_isCCW */
int GEOSCoordSeq_isCCW_r (
    GEOSContextHandle_t handle,
    const(GEOSCoordSequence)* s,
    char* is_ccw);

/* ========= Linear referencing functions ========= */

/** \see GEOSProject */
double GEOSProject_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* line,
    const(GEOSGeometry)* point);

/** \see GEOSInterpolate */
GEOSGeometry* GEOSInterpolate_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* line,
    double d);

/** \see GEOSProjectNormalized */
double GEOSProjectNormalized_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    const(GEOSGeometry)* p);

/** \see GEOSInterpolateNormalized */
GEOSGeometry* GEOSInterpolateNormalized_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double d);

/* ========== Buffer related functions ========== */

/** \see GEOSBuffer */
GEOSGeometry* GEOSBuffer_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double width,
    int quadsegs);

/**
* Cap styles control the ends of buffered lines.
* \see GEOSBuffer
*/
enum GEOSBufCapStyles
{
    /** End is rounded, with end point of original line in the centre of the round cap. */
    GEOSBUF_CAP_ROUND = 1,

    /** End is flat, with end point of original line at the end of the buffer */
    GEOSBUF_CAP_FLAT = 2,

    /** End is flat, with end point of original line in the middle of a square enclosing that point */
    GEOSBUF_CAP_SQUARE = 3
}

/**
* Join styles control the buffer shape at bends in a line.
* \see GEOSBuffer
*/
enum GEOSBufJoinStyles
{
    /**
    * Join is rounded, essentially each line is terminated
    * in a round cap. Form round corner.
    */
    GEOSBUF_JOIN_ROUND = 1,
    /**
    * Join is flat, with line between buffer edges,
    * through the join point. Forms flat corner.
    */
    GEOSBUF_JOIN_MITRE = 2,
    /**
    * Join is the point at which the two buffer edges intersect.
    * Forms sharp corner.
    */
    GEOSBUF_JOIN_BEVEL = 3
}

/** \see GEOSBufferParams_create */
GEOSBufferParams* GEOSBufferParams_create_r (GEOSContextHandle_t handle);

/** \see GEOSBufferParams_destroy */
void GEOSBufferParams_destroy_r (
    GEOSContextHandle_t handle,
    GEOSBufferParams* parms);

/** \see GEOSBufferParams_setEndCapStyle */
int GEOSBufferParams_setEndCapStyle_r (
    GEOSContextHandle_t handle,
    GEOSBufferParams* p,
    int style);

/** \see GEOSBufferParams_setJoinStyle */
int GEOSBufferParams_setJoinStyle_r (
    GEOSContextHandle_t handle,
    GEOSBufferParams* p,
    int joinStyle);

/** \see GEOSBufferParams_setMitreLimit */
int GEOSBufferParams_setMitreLimit_r (
    GEOSContextHandle_t handle,
    GEOSBufferParams* p,
    double mitreLimit);

/** \see GEOSBufferParams_setQuadrantSegments */
int GEOSBufferParams_setQuadrantSegments_r (
    GEOSContextHandle_t handle,
    GEOSBufferParams* p,
    int quadSegs);

/** \see GEOSBufferParams_setSingleSided */
int GEOSBufferParams_setSingleSided_r (
    GEOSContextHandle_t handle,
    GEOSBufferParams* p,
    int singleSided);

/** \see GEOSBufferWithParams */
GEOSGeometry* GEOSBufferWithParams_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    const(GEOSBufferParams)* p,
    double width);

/** \see GEOSBufferWithStyle */
GEOSGeometry* GEOSBufferWithStyle_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double width,
    int quadsegs,
    int endCapStyle,
    int joinStyle,
    double mitreLimit);

/** \see GEOSDensify */
GEOSGeometry* GEOSDensify_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double tolerance);

/** \see GEOSOffsetCurve */
GEOSGeometry* GEOSOffsetCurve_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double width,
    int quadsegs,
    int joinStyle,
    double mitreLimit);

/* ========= Geometry Constructors ========= */

/** \see GEOSGeom_createPoint */
GEOSGeometry* GEOSGeom_createPoint_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s);

/** \see GEOSGeom_createPointFromXY */
GEOSGeometry* GEOSGeom_createPointFromXY_r (
    GEOSContextHandle_t handle,
    double x,
    double y);

/** \see GEOSGeom_createEmptyPoint */
GEOSGeometry* GEOSGeom_createEmptyPoint_r (GEOSContextHandle_t handle);

/** \see GEOSGeom_createLinearRing */
GEOSGeometry* GEOSGeom_createLinearRing_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s);

/** \see GEOSGeom_createLineString */
GEOSGeometry* GEOSGeom_createLineString_r (
    GEOSContextHandle_t handle,
    GEOSCoordSequence* s);

/** \see GEOSGeom_createEmptyLineString */
GEOSGeometry* GEOSGeom_createEmptyLineString_r (GEOSContextHandle_t handle);

/** \see GEOSGeom_createEmptyPolygon */
GEOSGeometry* GEOSGeom_createEmptyPolygon_r (GEOSContextHandle_t handle);

/** \see GEOSGeom_createPolygon */
GEOSGeometry* GEOSGeom_createPolygon_r (
    GEOSContextHandle_t handle,
    GEOSGeometry* shell,
    GEOSGeometry** holes,
    uint nholes);

/** \see GEOSGeom_createCollection */
GEOSGeometry* GEOSGeom_createCollection_r (
    GEOSContextHandle_t handle,
    int type,
    GEOSGeometry** geoms,
    uint ngeoms);

/** \see GEOSGeom_createEmptyCollection */
GEOSGeometry* GEOSGeom_createEmptyCollection_r (
    GEOSContextHandle_t handle,
    int type);

/** \see GEOSGeom_clone */
GEOSGeometry* GEOSGeom_clone_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/* ========= Memory management ========= */

/** \see GEOSGeom_destroy */
void GEOSGeom_destroy_r (GEOSContextHandle_t handle, GEOSGeometry* g);

/* ========= Topology Operations ========= */

/** \see GEOSEnvelope */
GEOSGeometry* GEOSEnvelope_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSIntersection */
GEOSGeometry* GEOSIntersection_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSIntersectionPrec */
GEOSGeometry* GEOSIntersectionPrec_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double gridSize);

/** \see GEOSConvexHull */
GEOSGeometry* GEOSConvexHull_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSMinimumRotatedRectangle */
GEOSGeometry* GEOSMinimumRotatedRectangle_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSMaximumInscribedCircle */
GEOSGeometry* GEOSMaximumInscribedCircle_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double tolerance);

/** \see GEOSLargestEmptyCircle */
GEOSGeometry* GEOSLargestEmptyCircle_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    const(GEOSGeometry)* boundary,
    double tolerance);

/** \see GEOSMinimumWidth */
GEOSGeometry* GEOSMinimumWidth_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSMinimumClearanceLine */
GEOSGeometry* GEOSMinimumClearanceLine_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSMinimumClearance */
int GEOSMinimumClearance_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* distance);

/** \see GEOSDifference */
GEOSGeometry* GEOSDifference_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSDifferencePrec */
GEOSGeometry* GEOSDifferencePrec_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double gridSize);

/** \see GEOSSymDifference */
GEOSGeometry* GEOSSymDifference_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSSymDifferencePrec */
GEOSGeometry* GEOSSymDifferencePrec_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double gridSize);

/** \see GEOSBoundary */
GEOSGeometry* GEOSBoundary_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSUnion */
GEOSGeometry* GEOSUnion_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSUnionPrec */
GEOSGeometry* GEOSUnionPrec_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double gridSize);

/** \see GEOSUnaryUnion */
GEOSGeometry* GEOSUnaryUnion_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSUnaryUnionPrec */
GEOSGeometry* GEOSUnaryUnionPrec_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double gridSize);

/** \see GEOSCoverageUnion */
GEOSGeometry* GEOSCoverageUnion_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSPointOnSurface */
GEOSGeometry* GEOSPointOnSurface_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGetCentroid */
GEOSGeometry* GEOSGetCentroid_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSMinimumBoundingCircle */
GEOSGeometry* GEOSMinimumBoundingCircle_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* radius,
    GEOSGeometry** center);

/** \see GEOSNode */
GEOSGeometry* GEOSNode_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSClipByRect */
GEOSGeometry* GEOSClipByRect_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double xmin,
    double ymin,
    double xmax,
    double ymax);

/** \see GEOSPolygonize */
GEOSGeometry* GEOSPolygonize_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry*)* geoms,
    uint ngeoms);

/** \see GEOSPolygonize_valid */
GEOSGeometry* GEOSPolygonize_valid_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry*)* geoms,
    uint ngems);

/** \see GEOSPolygonizer_getCutEdges */
GEOSGeometry* GEOSPolygonizer_getCutEdges_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry*)* geoms,
    uint ngeoms);

/** \see GEOSPolygonize_full */
GEOSGeometry* GEOSPolygonize_full_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* input,
    GEOSGeometry** cuts,
    GEOSGeometry** dangles,
    GEOSGeometry** invalidRings);

/** \see GEOSBuildArea */
GEOSGeometry* GEOSBuildArea_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSLineMerge */
GEOSGeometry* GEOSLineMerge_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSReverse */
GEOSGeometry* GEOSReverse_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSSimplify */
GEOSGeometry* GEOSSimplify_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double tolerance);

/** \see GEOSTopologyPreserveSimplify */
GEOSGeometry* GEOSTopologyPreserveSimplify_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double tolerance);

/** \see GEOSGeom_extractUniquePoints */
GEOSGeometry* GEOSGeom_extractUniquePoints_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSSharedPaths */
GEOSGeometry* GEOSSharedPaths_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSSnap */
GEOSGeometry* GEOSSnap_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double tolerance);

/** \see GEOSDelaunayTriangulation */
GEOSGeometry* GEOSDelaunayTriangulation_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double tolerance,
    int onlyEdges);

/** \see GEOSConstrainedDelaunayTriangulation */
GEOSGeometry* GEOSConstrainedDelaunayTriangulation_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSVoronoiDiagram */
GEOSGeometry* GEOSVoronoiDiagram_r (
    GEOSContextHandle_t extHandle,
    const(GEOSGeometry)* g,
    const(GEOSGeometry)* env,
    double tolerance,
    int onlyEdges);

/** \see GEOSSegmentIntersection */
int GEOSSegmentIntersection_r (
    GEOSContextHandle_t extHandle,
    double ax0,
    double ay0,
    double ax1,
    double ay1,
    double bx0,
    double by0,
    double bx1,
    double by1,
    double* cx,
    double* cy);

/* ========= Binary predicates ========= */

/** \see GEOSDisjoint */
char GEOSDisjoint_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSTouches */
char GEOSTouches_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSIntersects */
char GEOSIntersects_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSCrosses */
char GEOSCrosses_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSWithin */
char GEOSWithin_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSContains */
char GEOSContains_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSOverlaps */
char GEOSOverlaps_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSEquals */
char GEOSEquals_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSEqualsExact */
char GEOSEqualsExact_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double tolerance);

/** \see GEOSCovers */
char GEOSCovers_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSCoveredBy */
char GEOSCoveredBy_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/* ========= Prepared Geometry Binary Predicates ========== */

/** \see GEOSPrepare */
const(GEOSPreparedGeometry)* GEOSPrepare_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSPreparedGeom_destroy */
void GEOSPreparedGeom_destroy_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* g);

/** \see GEOSPreparedContains */
char GEOSPreparedContains_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedContainsProperly */
char GEOSPreparedContainsProperly_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedCoveredBy */
char GEOSPreparedCoveredBy_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedCovers */
char GEOSPreparedCovers_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedCrosses */
char GEOSPreparedCrosses_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedDisjoint */
char GEOSPreparedDisjoint_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedIntersects */
char GEOSPreparedIntersects_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedOverlaps */
char GEOSPreparedOverlaps_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedTouches */
char GEOSPreparedTouches_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedWithin */
char GEOSPreparedWithin_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedNearestPoints */
GEOSCoordSequence* GEOSPreparedNearestPoints_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/** \see GEOSPreparedDistance */
int GEOSPreparedDistance_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2,
    double* dist);

/** \see GEOSPreparedDistanceWithin */
char GEOSPreparedDistanceWithin_r (
    GEOSContextHandle_t handle,
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2,
    double dist);

/* ========== STRtree ========== */

/** \see GEOSSTRtree_create */
GEOSSTRtree* GEOSSTRtree_create_r (
    GEOSContextHandle_t handle,
    size_t nodeCapacity);

/** \see GEOSSTRtree_insert */
void GEOSSTRtree_insert_r (
    GEOSContextHandle_t handle,
    GEOSSTRtree* tree,
    const(GEOSGeometry)* g,
    void* item);

/** \see GEOSSTRtree_query */
void GEOSSTRtree_query_r (
    GEOSContextHandle_t handle,
    GEOSSTRtree* tree,
    const(GEOSGeometry)* g,
    GEOSQueryCallback callback,
    void* userdata);

/** \see GEOSSTRtree_nearest */
const(GEOSGeometry)* GEOSSTRtree_nearest_r (
    GEOSContextHandle_t handle,
    GEOSSTRtree* tree,
    const(GEOSGeometry)* geom);

/** \see GEOSSTRtree_nearest_generic */
const(void)* GEOSSTRtree_nearest_generic_r (
    GEOSContextHandle_t handle,
    GEOSSTRtree* tree,
    const(void)* item,
    const(GEOSGeometry)* itemEnvelope,
    GEOSDistanceCallback distancefn,
    void* userdata);

/** \see GEOSSTRtree_iterate */
void GEOSSTRtree_iterate_r (
    GEOSContextHandle_t handle,
    GEOSSTRtree* tree,
    GEOSQueryCallback callback,
    void* userdata);

/** \see GEOSSTRtree_remove */
char GEOSSTRtree_remove_r (
    GEOSContextHandle_t handle,
    GEOSSTRtree* tree,
    const(GEOSGeometry)* g,
    void* item);

/** \see GEOSSTRtree_destroy */
void GEOSSTRtree_destroy_r (GEOSContextHandle_t handle, GEOSSTRtree* tree);

/* ========= Unary predicate ========= */

/** \see GEOSisEmpty */
char GEOSisEmpty_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSisSimple */
char GEOSisSimple_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSisRing */
char GEOSisRing_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSHasZ */
char GEOSHasZ_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSisClosed */
char GEOSisClosed_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/* ========== Dimensionally Extended 9 Intersection Model ========== */

/**
* Controls the behavior of the result of GEOSRelate when returning
* DE9IM results for two geometries.
*/
enum GEOSRelateBoundaryNodeRules
{
    /** See geos::algorithm::BoundaryNodeRule::getBoundaryRuleMod2() */
    GEOSRELATE_BNR_MOD2 = 1,
    /** Same as \ref GEOSRELATE_BNR_MOD2 */
    GEOSRELATE_BNR_OGC = 1,
    /** See geos::algorithm::BoundaryNodeRule::getBoundaryEndPoint() */
    GEOSRELATE_BNR_ENDPOINT = 2,
    /** See geos::algorithm::BoundaryNodeRule::getBoundaryMultivalentEndPoint() */
    GEOSRELATE_BNR_MULTIVALENT_ENDPOINT = 3,
    /** See geos::algorithm::BoundaryNodeRule::getBoundaryMonovalentEndPoint() */
    GEOSRELATE_BNR_MONOVALENT_ENDPOINT = 4
}

/** \see GEOSRelatePattern */
char GEOSRelatePattern_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    const(char)* pat);

/** \see GEOSRelate */
char* GEOSRelate_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/** \see GEOSRelatePatternMatch */
char GEOSRelatePatternMatch_r (
    GEOSContextHandle_t handle,
    const(char)* mat,
    const(char)* pat);

/** \see GEOSRelateBoundaryNodeRule */
char* GEOSRelateBoundaryNodeRule_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    int bnr);

/* ========= Validity checking ========= */

/** Change behaviour of validity testing in \ref GEOSisValidDetail */
enum GEOSValidFlags
{
    /** Allow self-touching rings to form a hole in a polygon. */
    GEOSVALID_ALLOW_SELFTOUCHING_RING_FORMING_HOLE = 1
}

/** \see GEOSisValid */
char GEOSisValid_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSisValidReason */
char* GEOSisValidReason_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSisValidDetail */
char GEOSisValidDetail_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    int flags,
    char** reason,
    GEOSGeometry** location);

/* ========== Make Valid ========== */

/**
* Algorithm to use when repairing invalid geometries.
*
* \see GEOSMakeValidWithParams
*/
enum GEOSMakeValidMethods
{
    /** Original method, combines all rings into
        a set of noded lines and then extracts valid
        polygons from that linework. */
    GEOS_MAKE_VALID_LINEWORK = 0,
    /** Structured method, first makes all rings valid
        then merges shells and subtracts holes from
        shells to generate valid result. Assumes that
        holes and shells are correctly categorized. */
    GEOS_MAKE_VALID_STRUCTURE = 1
}

/** \see GEOSMakeValidParams_create */
GEOSMakeValidParams* GEOSMakeValidParams_create_r (
    GEOSContextHandle_t extHandle);

/** \see GEOSMakeValidParams_destroy */
void GEOSMakeValidParams_destroy_r (
    GEOSContextHandle_t handle,
    GEOSMakeValidParams* parms);

/** \see GEOSMakeValidParams_setKeepCollapsed */
int GEOSMakeValidParams_setKeepCollapsed_r (
    GEOSContextHandle_t handle,
    GEOSMakeValidParams* p,
    int style);

/** \see GEOSMakeValidParams_setMethod */
int GEOSMakeValidParams_setMethod_r (
    GEOSContextHandle_t handle,
    GEOSMakeValidParams* p,
    GEOSMakeValidMethods method);

/** \see GEOSMakeValid */
GEOSGeometry* GEOSMakeValid_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSMakeValidWithParams */
GEOSGeometry* GEOSMakeValidWithParams_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    const(GEOSMakeValidParams)* makeValidParams);

/* ========== Geometry info ========== */

/** \see GEOSGeomType */
/* Return NULL on exception, result must be freed by caller. */
char* GEOSGeomType_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSGeomTypeId */
int GEOSGeomTypeId_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSGetSRID */
int GEOSGetSRID_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSSetSRID */
void GEOSSetSRID_r (GEOSContextHandle_t handle, GEOSGeometry* g, int SRID);

/** \see GEOSGeom_getUserData */
void* GEOSGeom_getUserData_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGeom_setUserData */
void GEOSGeom_setUserData_r (
    GEOSContextHandle_t handle,
    GEOSGeometry* g,
    void* userData);

/** \see GEOSGetNumGeometries */
int GEOSGetNumGeometries_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSGetGeometryN */
const(GEOSGeometry)* GEOSGetGeometryN_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    int n);

/** \see GEOSNormalize */
int GEOSNormalize_r (GEOSContextHandle_t handle, GEOSGeometry* g);

/**
* Controls the behavior of GEOSGeom_setPrecision()
* when altering the precision of a geometry.
*/
enum GEOSPrecisionRules
{
    /** The output is always valid. Collapsed geometry elements (including both polygons and lines) are removed. */
    GEOS_PREC_VALID_OUTPUT = 0,
    /** Precision reduction is performed pointwise. Output geometry may be invalid due to collapse or self-intersection. (This might be better called "GEOS_PREC_POINTWISE" - the current name is historical.) */
    GEOS_PREC_NO_TOPO = 1,
    /** Like the default mode, except that collapsed linear geometry elements are preserved. Collapsed polygonal input elements are removed. */
    GEOS_PREC_KEEP_COLLAPSED = 2
}

/** \see GEOSGeom_setPrecision */
GEOSGeometry* GEOSGeom_setPrecision_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double gridSize,
    int flags);

/** \see GEOSGeom_getPrecision */
double GEOSGeom_getPrecision_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGetNumInteriorRings */
int GEOSGetNumInteriorRings_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGeomGetNumPoints */
int GEOSGeomGetNumPoints_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/** \see GEOSGeomGetX */
int GEOSGeomGetX_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* x);

/** \see GEOSGeomGetY */
int GEOSGeomGetY_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* y);

/** \see GEOSGeomGetZ */
int GEOSGeomGetZ_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* z);

/** \see GEOSGetInteriorRingN */
const(GEOSGeometry)* GEOSGetInteriorRingN_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    int n);

/** \see GEOSGetExteriorRing */
const(GEOSGeometry)* GEOSGetExteriorRing_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGetNumCoordinates */
int GEOSGetNumCoordinates_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGeom_getCoordSeq */
const(GEOSCoordSequence)* GEOSGeom_getCoordSeq_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGeom_getDimensions */
int GEOSGeom_getDimensions_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGeom_getCoordinateDimension */
int GEOSGeom_getCoordinateDimension_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGeom_getXMin */
int GEOSGeom_getXMin_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* value);

/** \see GEOSGeom_getYMin */
int GEOSGeom_getYMin_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* value);

/** \see GEOSGeom_getXMax */
int GEOSGeom_getXMax_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* value);

/** \see GEOSGeom_getYMax */
int GEOSGeom_getYMax_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* value);

/** \see GEOSGeomGetPointN */
GEOSGeometry* GEOSGeomGetPointN_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    int n);

/** \see GEOSGeomGetStartPoint */
GEOSGeometry* GEOSGeomGetStartPoint_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/** \see GEOSGeomGetEndPoint */
GEOSGeometry* GEOSGeomGetEndPoint_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/* ========= Misc functions ========= */

/** \see GEOSArea */
int GEOSArea_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* area);

/** \see GEOSLength */
int GEOSLength_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* length);

/** \see GEOSDistance */
int GEOSDistance_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/** \see GEOSDistanceWithin */
char GEOSDistanceWithin_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double dist);

/** \see GEOSDistanceIndexed */
int GEOSDistanceIndexed_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/** \see GEOSHausdorffDistance */
int GEOSHausdorffDistance_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/** \see GEOSHausdorffDistanceDensify */
int GEOSHausdorffDistanceDensify_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double densifyFrac,
    double* dist);

/** \see GEOSFrechetDistance */
int GEOSFrechetDistance_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/** \see GEOSFrechetDistanceDensify */
int GEOSFrechetDistanceDensify_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double densifyFrac,
    double* dist);

/** \see GEOSGeomGetLength */
int GEOSGeomGetLength_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double* length);

/** \see GEOSNearestPoints */
GEOSCoordSequence* GEOSNearestPoints_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/* ========= Algorithms ========= */

/** \see GEOSOrientationIndex */
int GEOSOrientationIndex_r (
    GEOSContextHandle_t handle,
    double Ax,
    double Ay,
    double Bx,
    double By,
    double Px,
    double Py);

/* ========== Reader and Writer APIs ========== */

/**
* Reader object to read Well-Known Text (WKT) format and construct Geometry.
* \see GEOSWKTReader_create
* \see GEOSWKTReader_create_r
*/
struct GEOSWKTReader_t;
alias GEOSWKTReader = GEOSWKTReader_t;

/**
* Writer object to turn Geometry into Well-Known Text (WKT).
* \see GEOSWKTWriter_create
* \see GEOSWKTWriter_create_r
*/
struct GEOSWKTWriter_t;
alias GEOSWKTWriter = GEOSWKTWriter_t;

/**
* Reader object to read Well-Known Binary (WKB) format and construct Geometry.
* \see GEOSWKBReader_create
* \see GEOSWKBReader_create_r
*/
struct GEOSWKBReader_t;
alias GEOSWKBReader = GEOSWKBReader_t;

/**
* Writer object to turn Geometry into Well-Known Binary (WKB).
* \see GEOSWKBWriter_create
* \see GEOSWKBWriter_create_r
*/
struct GEOSWKBWriter_t;
alias GEOSWKBWriter = GEOSWKBWriter_t;

/**
* Reader object to read GeoJSON format and construct a Geometry.
* \see GEOSGeoJSONReader_create
* \see GEOSGeoJSONReader_create_r
*/
struct GEOSGeoJSONReader_t;
alias GEOSGeoJSONReader = GEOSGeoJSONReader_t;

/**
* Writer object to turn a Geometry into GeoJSON.
* \see GEOSGeoJSONReader_create
* \see GEOSGeoJSONReader_create_r
*/
struct GEOSGeoJSONWriter_t;
alias GEOSGeoJSONWriter = GEOSGeoJSONWriter_t;

/* ========== WKT Reader ========== */

/** \see GEOSWKTReader_create */
GEOSWKTReader* GEOSWKTReader_create_r (GEOSContextHandle_t handle);

/** \see GEOSWKTReader_destroy */
void GEOSWKTReader_destroy_r (
    GEOSContextHandle_t handle,
    GEOSWKTReader* reader);

/** \see GEOSWKTReader_read */
GEOSGeometry* GEOSWKTReader_read_r (
    GEOSContextHandle_t handle,
    GEOSWKTReader* reader,
    const(char)* wkt);

/* ========== WKT Writer ========== */

/** \see GEOSWKTReader_create */
GEOSWKTWriter* GEOSWKTWriter_create_r (GEOSContextHandle_t handle);

/** \see GEOSWKTWriter_destroy */
void GEOSWKTWriter_destroy_r (
    GEOSContextHandle_t handle,
    GEOSWKTWriter* writer);

/** \see GEOSWKTWriter_write */
char* GEOSWKTWriter_write_r (
    GEOSContextHandle_t handle,
    GEOSWKTWriter* writer,
    const(GEOSGeometry)* g);

/** \see GEOSWKTWriter_setTrim */
void GEOSWKTWriter_setTrim_r (
    GEOSContextHandle_t handle,
    GEOSWKTWriter* writer,
    char trim);

/** \see GEOSWKTWriter_setRoundingPrecision */
void GEOSWKTWriter_setRoundingPrecision_r (
    GEOSContextHandle_t handle,
    GEOSWKTWriter* writer,
    int precision);

/** \see GEOSWKTWriter_setOutputDimension */
void GEOSWKTWriter_setOutputDimension_r (
    GEOSContextHandle_t handle,
    GEOSWKTWriter* writer,
    int dim);

/** \see GEOSWKTWriter_getOutputDimension */
int GEOSWKTWriter_getOutputDimension_r (
    GEOSContextHandle_t handle,
    GEOSWKTWriter* writer);

/** \see GEOSWKTWriter_setOld3D */
void GEOSWKTWriter_setOld3D_r (
    GEOSContextHandle_t handle,
    GEOSWKTWriter* writer,
    int useOld3D);

/* ========== WKB Reader ========== */

/** \see GEOSWKBReader_create */
GEOSWKBReader* GEOSWKBReader_create_r (GEOSContextHandle_t handle);

/** \see GEOSWKBReader_destroy */
void GEOSWKBReader_destroy_r (
    GEOSContextHandle_t handle,
    GEOSWKBReader* reader);

/** \see GEOSWKBReader_read */
GEOSGeometry* GEOSWKBReader_read_r (
    GEOSContextHandle_t handle,
    GEOSWKBReader* reader,
    const(ubyte)* wkb,
    size_t size);

/** \see GEOSWKBReader_readHEX */
GEOSGeometry* GEOSWKBReader_readHEX_r (
    GEOSContextHandle_t handle,
    GEOSWKBReader* reader,
    const(ubyte)* hex,
    size_t size);

/* ========== WKB Writer ========== */

/** \see GEOSWKBWriter_create */
GEOSWKBWriter* GEOSWKBWriter_create_r (GEOSContextHandle_t handle);

/** \see GEOSWKBWriter_destroy */
void GEOSWKBWriter_destroy_r (
    GEOSContextHandle_t handle,
    GEOSWKBWriter* writer);

/** \see GEOSWKBWriter_write */
ubyte* GEOSWKBWriter_write_r (
    GEOSContextHandle_t handle,
    GEOSWKBWriter* writer,
    const(GEOSGeometry)* g,
    size_t* size);

/** \see GEOSWKBWriter_writeHEX */
ubyte* GEOSWKBWriter_writeHEX_r (
    GEOSContextHandle_t handle,
    GEOSWKBWriter* writer,
    const(GEOSGeometry)* g,
    size_t* size);

/** \see GEOSWKBWriter_getOutputDimension */
int GEOSWKBWriter_getOutputDimension_r (
    GEOSContextHandle_t handle,
    const(GEOSWKBWriter)* writer);

/** \see GEOSWKBWriter_setOutputDimension */
void GEOSWKBWriter_setOutputDimension_r (
    GEOSContextHandle_t handle,
    GEOSWKBWriter* writer,
    int newDimension);

/** \see GEOSWKBWriter_getByteOrder */
int GEOSWKBWriter_getByteOrder_r (
    GEOSContextHandle_t handle,
    const(GEOSWKBWriter)* writer);

/** \see GEOSWKBWriter_setByteOrder */
void GEOSWKBWriter_setByteOrder_r (
    GEOSContextHandle_t handle,
    GEOSWKBWriter* writer,
    int byteOrder);

/** \see GEOSWKBWriter_getFlavor */
int GEOSWKBWriter_getFlavor_r (
    GEOSContextHandle_t handle,
    const(GEOSWKBWriter)* writer);

/** \see GEOSWKBWriter_setFlavor */
void GEOSWKBWriter_setFlavor_r (
    GEOSContextHandle_t handle,
    GEOSWKBWriter* writer,
    int flavor);

/** \see GEOSWKBWriter_getIncludeSRID */
char GEOSWKBWriter_getIncludeSRID_r (
    GEOSContextHandle_t handle,
    const(GEOSWKBWriter)* writer);

/** \see GEOSWKBWriter_setIncludeSRID */
void GEOSWKBWriter_setIncludeSRID_r (
    GEOSContextHandle_t handle,
    GEOSWKBWriter* writer,
    const char writeSRID);

/* ========== GeoJSON Reader ========== */

/** \see GEOSGeoJSONReader_create */
GEOSGeoJSONReader* GEOSGeoJSONReader_create_r (GEOSContextHandle_t handle);

/** \see GEOSGeoJSONReader_destroy */
void GEOSGeoJSONReader_destroy_r (
    GEOSContextHandle_t handle,
    GEOSGeoJSONReader* reader);

/** \see GEOSWKTReader_read */
GEOSGeometry* GEOSGeoJSONReader_readGeometry_r (
    GEOSContextHandle_t handle,
    GEOSGeoJSONReader* reader,
    const(char)* geojson);

/* ========== GeoJSON Writer ========== */

/** \see GEOSGeoJSONWriter_create */
GEOSGeoJSONWriter* GEOSGeoJSONWriter_create_r (GEOSContextHandle_t handle);

/** \see GEOSGeoJSONWriter_destroy */
void GEOSGeoJSONWriter_destroy_r (
    GEOSContextHandle_t handle,
    GEOSGeoJSONWriter* writer);

/** \see GEOSGeoJSONWriter_writeGeometry */
char* GEOSGeoJSONWriter_writeGeometry_r (
    GEOSContextHandle_t handle,
    GEOSGeoJSONWriter* writer,
    const(GEOSGeometry)* g,
    int indent);

/** \see GEOSFree */
void GEOSFree_r (GEOSContextHandle_t handle, void* buffer);

/**
* Returns the current GEOS version string. eg: "3.10.0"
* This function does not have a reentrant variant and is
* available if `GEOS_USE_ONLY_R_API` is defined.
* \return version string
*/
const(char)* GEOSversion ();

/*
* External code to GEOS can define GEOS_USE_ONLY_R_API
* to strip the non-reentrant API functions from this header,
* leaving only the "_r" compatible variants.
*/

/* ========== Initialization, cleanup, version ========== */

/**
* For non-reentrant code, set up an execution contact, and associate
* \ref GEOSMessageHandler functions with it, to pass error and notice
* messages back to the calling application.
* <pre>
* typedef void (*GEOSMessageHandler)(const char *fmt, ...);
* </pre>
*
* \param notice_function Handle notice messages
* \param error_function Handle error messages
*/
void initGEOS (
    GEOSMessageHandler notice_function,
    GEOSMessageHandler error_function);

/**
* For non-reentrant code, call when all GEOS operations are complete,
* cleans up global resources.
*/
void finishGEOS ();

/* ========= Coordinate Sequence functions ========= */

/**
* Create a coordinate sequence.
* \param size number of coordinates in the sequence
* \param dims dimensionality of the coordinates (2 or 3)
* \return the sequence or NULL on exception
*/
GEOSCoordSequence* GEOSCoordSeq_create (uint size, uint dims);

/**
* Create a coordinate sequence by copying from a buffer of doubles (XYXY or XYZXYZ)
* \param buf pointer to buffer
* \param size number of coordinates in the sequence
* \param hasZ does buffer have Z values?
* \param hasM does buffer have M values? (they will be ignored)
* \return the sequence or NULL on exception
*/
GEOSCoordSequence* GEOSCoordSeq_copyFromBuffer (const(double)* buf, uint size, int hasZ, int hasM);

/**
* Create a coordinate sequence by copying from arrays of doubles
* \param x array of x coordinates
* \param y array of y coordinates
* \param z array of z coordinates, or NULL
* \param m array of m coordinates, (must be NULL)
* \param size length of each array
* \return the sequence or NULL on exception
*/
GEOSCoordSequence* GEOSCoordSeq_copyFromArrays (const(double)* x, const(double)* y, const(double)* z, const(double)* m, uint size);

/**
* Copy the contents of a coordinate sequence to a buffer of doubles (XYXY or XYZXYZ)
* \param s sequence to copy
* \param buf buffer to which coordinates should be copied
* \param hasZ copy Z values to buffer?
* \param hasM copy M values to buffer? (will be NaN)
* \return 1 on success, 0 on error
*/
int GEOSCoordSeq_copyToBuffer (const(GEOSCoordSequence)* s, double* buf, int hasZ, int hasM);

/**
* Copy the contents of a coordinate sequence to a buffer of doubles (XYZY or XYZXYZ)
* \param s sequence to copy
* \param x array to which x values should be copied
* \param y array to which y values should be copied
* \param z array to which z values should be copied, or NULL
* \param m array to which m values should be copied (will all be NAN)
* \return 1 on success, 0 on error
*/
int GEOSCoordSeq_copyToArrays (const(GEOSCoordSequence)* s, double* x, double* y, double* z, double* m);

/**
* Clone a coordinate sequence.
* \param s the coordinate sequence to clone
* \return a copy of the coordinate sequence or NULL on exception
*/
GEOSCoordSequence* GEOSCoordSeq_clone (const(GEOSCoordSequence)* s);

/**
* Destroy a coordinate sequence, freeing all memory.
* \param s the coordinate sequence to destroy
*/
void GEOSCoordSeq_destroy (GEOSCoordSequence* s);

/**
* Set X ordinate values in a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param val the value to set the ordinate to
* \return 0 on exception
*/
int GEOSCoordSeq_setX (GEOSCoordSequence* s, uint idx, double val);
/**
* Set Y ordinate values in a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param val the value to set the ordinate to
* \return 0 on exception
*/
int GEOSCoordSeq_setY (GEOSCoordSequence* s, uint idx, double val);
/**
* Set Z ordinate values in a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param val the value to set the ordinate to
* \return 0 on exception
*/
int GEOSCoordSeq_setZ (GEOSCoordSequence* s, uint idx, double val);
/**
* Set X and Y ordinate values in a coordinate sequence simultaneously.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param x the value to set the X ordinate to
* \param y the value to set the Y ordinate to
* \return 0 on exception
*/
int GEOSCoordSeq_setXY (GEOSCoordSequence* s, uint idx, double x, double y);
/**
* Set X, Y and Z ordinate values in a coordinate sequence simultaneously.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param x the value to set the X ordinate to
* \param y the value to set the Y ordinate to
* \param z the value to set the Z ordinate to
* \return 0 on exception
*/
int GEOSCoordSeq_setXYZ (
    GEOSCoordSequence* s,
    uint idx,
    double x,
    double y,
    double z);
/**
* Set Nth ordinate value in a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param dim the dimension number of the ordinate to alter, zero based
* \param val the value to set the ordinate to
* \return 0 on exception
*/
int GEOSCoordSeq_setOrdinate (
    GEOSCoordSequence* s,
    uint idx,
    uint dim,
    double val);

/**
* Read X ordinate values from a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param val pointer where ordinate value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getX (const(GEOSCoordSequence)* s, uint idx, double* val);

/**
* Read Y ordinate values from a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param val pointer where ordinate value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getY (const(GEOSCoordSequence)* s, uint idx, double* val);
/**
* Read Z ordinate values from a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param val pointer where ordinate value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getZ (const(GEOSCoordSequence)* s, uint idx, double* val);
/**
* Read X and Y ordinate values from a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param x pointer where ordinate X value will be placed
* \param y pointer where ordinate Y value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getXY (
    const(GEOSCoordSequence)* s,
    uint idx,
    double* x,
    double* y);
/**
* Read X and Y ordinate values from a coordinate sequence.
* \param s the coordinate sequence
* \param idx the index of the coordinate to alter, zero based
* \param x pointer where ordinate X value will be placed
* \param y pointer where ordinate Y value will be placed
* \param z pointer where ordinate Z value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getXYZ (
    const(GEOSCoordSequence)* s,
    uint idx,
    double* x,
    double* y,
    double* z);
/**
* Read Nth ordinate value from a coordinate sequence.
* \param[in] s the coordinate sequence
* \param[in] idx the index of the coordinate to alter, zero based
* \param[in] dim the dimension number of the ordinate to read, zero based
* \param[out] val pointer where ordinate value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getOrdinate (
    const(GEOSCoordSequence)* s,
    uint idx,
    uint dim,
    double* val);

/**
* Get size info from a coordinate sequence.
* \param[in] s the coordinate sequence
* \param[out] size pointer where size value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getSize (const(GEOSCoordSequence)* s, uint* size);

/**
* Get dimension info from a coordinate sequence.
* \param[in] s the coordinate sequence
* \param[out] dims pointer where dimension value will be placed
* \return 0 on exception
*/
int GEOSCoordSeq_getDimensions (const(GEOSCoordSequence)* s, uint* dims);

/**
* Check orientation of a coordinate sequence. Closure of the sequence is
* assumed. Invalid (collapsed) sequences will return false. Short (less
* than 4 points) sequences will return exception.
* \param s the coordinate sequence
* \param is_ccw pointer for ccw value, 1 if counter-clockwise orientation, 0 otherwise
* \return 0 on exception, 1 on success
*/
int GEOSCoordSeq_isCCW (const(GEOSCoordSequence)* s, char* is_ccw);

/* ========== Linear referencing functions */

/**
* Distance of point projected onto line from the start of the line.
* \param line linear target of projection
* \param point point to be projected onto 'g'
* \return distance along line that point projects to, -1 on exception
*
* \note Line parameter must be a LineString.
*/
double GEOSProject (const(GEOSGeometry)* line, const(GEOSGeometry)* point);

/**
* Measuring from start of line, return point that is distance
* the start. Line parameter must be a LineString.
* \param line linear target of projection
* \param d distance from start of line to created point
* \return The point \ref GEOSGeometry that is distance from the start of line.
* Caller takes ownership of returned geometry.
*/
GEOSGeometry* GEOSInterpolate (const(GEOSGeometry)* line, double d);

/**
* Project point to line and calculate the **proportion** of
* the line the point is from the start. For example, a point that
* projects to the middle of a line would be return 0.5.
* \param line linear target of projection
* \param point the point to project
* \return The proportion of the overall line length that the projected
* point falls at.
*/
double GEOSProjectNormalized (
    const(GEOSGeometry)* line,
    const(GEOSGeometry)* point);

/**
* Measuring from start of line, return point that is a proportion
* the start. Line parameter must be a LineString.
* \param line linear target of projection
* \param proportion The proportion from the start of line to created point
* \return The point \ref GEOSGeometry that is distance from the start of line.
* Caller takes ownership of returned geometry.
*/
GEOSGeometry* GEOSInterpolateNormalized (
    const(GEOSGeometry)* line,
    double proportion);

/* ========== Buffer related functions ========== */

/**
* Buffer a geometry.
* \param g The input geometry to be buffered.
* \param width The distance by which to expand the geometry (or contract)
*        if the value is negative.
* \param quadsegs The number of segments per quadrant to generate. More
*        segments provides a more "precise" buffer at the expense of size.
* \return A \ref GEOSGeometry of the buffered result.
* NULL on exception. Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSBuffer (const(GEOSGeometry)* g, double width, int quadsegs);

/**
* Create a default GEOSBufferParams object for controlling the shape
* of buffered generated by \ref GEOSBuffer.
* \return A newly allocated GEOSBufferParams. NULL on exception.
* Caller is responsible for freeing with GEOSBufferParams_destroy().
*/
GEOSBufferParams* GEOSBufferParams_create ();

/**
* Destroy a GEOSBufferParams and free all associated memory.
* \param parms The object to destroy.
*/
void GEOSBufferParams_destroy (GEOSBufferParams* parms);

/**
* Set the end cap type of a GEOSBufferParams to the desired style,
* which must be one enumerated in \ref GEOSBufCapStyles.
* \return 0 on exception, 1 on success.
*/
int GEOSBufferParams_setEndCapStyle (GEOSBufferParams* p, int style);

/**
* Set the join type of a GEOSBufferParams to the desired style,
* which must be one enumerated in \ref GEOSBufJoinStyles.
* \return 0 on exception, 1 on success.
*/
int GEOSBufferParams_setJoinStyle (GEOSBufferParams* p, int joinStyle);

/**
* Set the mitre limit of a GEOSBufferParams to the desired size.
* For acute angles, a mitre join can extend very very far from
* the input geometry, which is probably not desired. The
* mitre limit places an upper bound on that.
* \param p The GEOSBufferParams to operate on
* \param mitreLimit The limit to set
* \return 0 on exception, 1 on success.
*/
int GEOSBufferParams_setMitreLimit (GEOSBufferParams* p, double mitreLimit);

/**
* Set the number of segments to use to stroke each quadrant
* of circular arcs generated by the buffering process. More
* segments means a smoother output, but with larger size.
* \param p The GEOSBufferParams to operate on
* \param quadSegs Number of segments per quadrant
* \return 0 on exception, 1 on success.
*/
int GEOSBufferParams_setQuadrantSegments (GEOSBufferParams* p, int quadSegs);

/**
* Sets whether the computed buffer should be single-sided.
* A single-sided buffer is constructed on only one side of each input line.
* \see geos::operation::buffer::BufferParameters::setSingleSided
* \param p The GEOSBufferParams to operate on
* \param singleSided Set to 1 for single-sided output 0 otherwise
* \return 0 on exception, 1 on success.
*/
int GEOSBufferParams_setSingleSided (GEOSBufferParams* p, int singleSided);

/**
* Generates a buffer using the special parameters in the GEOSBufferParams
* \param g The geometry to buffer
* \param p The parameters to apply to the buffer process
* \param width The buffer distance
* \return The buffered geometry, or NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSBufferWithParams (
    const(GEOSGeometry)* g,
    const(GEOSBufferParams)* p,
    double width);

/**
* Generate a buffer using the provided style parameters.
* \param g The geometry to buffer
* \param width Width of the buffer
* \param quadsegs Number of segments per quadrant
* \param endCapStyle See \ref GEOSBufCapStyles
* \param joinStyle See \ref GEOSBufJoinStyles
* \param mitreLimit See GEOSBufferParams_setMitreLimit
* \return The buffered geometry, or NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSBufferWithStyle (
    const(GEOSGeometry)* g,
    double width,
    int quadsegs,
    int endCapStyle,
    int joinStyle,
    double mitreLimit);

/**
* Densifies a geometry using a given distance tolerance.
* Additional vertices will be added to every line segment
* that is greater this tolerance; these vertices will
* evenly subdivide that segment.
* Only linear components of input geometry are densified.
* \param g The geometry to densify
* \param tolerance the distance tolerance to densify
* \return The densified geometry, or NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSDensify (const(GEOSGeometry)* g, double tolerance);

/**
* Generates offset curve for linear geometry.
* Only LineStrings accepted as input.
* \param g The linear geometry to offset from
* \param width Distance to offset from the curve.
*        Negative for a right-side offset.
*        Positive for a left-side offset.
* \param quadsegs Number of segments per quadrant
* \param joinStyle See \ref GEOSBufJoinStyles
* \param mitreLimit See GEOSBufferParams_setMitreLimit
* \return The offset geometry. Returns NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::buffer::BufferBuilder::bufferLineSingleSided
*/
GEOSGeometry* GEOSOffsetCurve (
    const(GEOSGeometry)* g,
    double width,
    int quadsegs,
    int joinStyle,
    double mitreLimit);

/* ========= Geometry Constructors ========= */

/**
* Creates a point geometry from a coordinate sequence.
* \param s Input coordinate sequence, ownership passes to the geometry
* \return A newly allocated point geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createPoint (GEOSCoordSequence* s);

/**
* Creates a point geometry from a pair of coordinates.
* \param x The X coordinate
* \param y The Y coordinate
* \return A newly allocated point geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createPointFromXY (double x, double y);

/**
* Creates an empty point.
* \return A newly allocated empty point geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createEmptyPoint ();

/**
* Creates a linear ring geometry, for use in a polygon.
* \param s Input coordinate sequence, ownership passes to the geometry
* \return A newly allocated linear ring geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createLinearRing (GEOSCoordSequence* s);

/**
* Creates a linestring geometry.
* \param s Input coordinate sequence, ownership passes to the geometry
* \return A newly allocated linestring geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createLineString (GEOSCoordSequence* s);

/**
* Creates an emptylinestring geometry.
* \return A newly allocated linestring geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createEmptyLineString ();

/**
* Creates an empty polygon geometry.
* \return A newly allocated empty polygon geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createEmptyPolygon ();

/**
* Creates a polygon geometry from line ring geometries.
* \param shell A linear ring that is the exterior ring of the polygon.
* \param holes An array of linear rings that are the holes.
* \param nholes The number of rings in the holes array.
* \return A newly allocated geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \note The holes argument is an array of GEOSGeometry* objects.
*       The caller **retains ownership** of the containing array,
*       but the ownership of the pointed-to objects is transferred
*       to the returned \ref GEOSGeometry.
*/
GEOSGeometry* GEOSGeom_createPolygon (
    GEOSGeometry* shell,
    GEOSGeometry** holes,
    uint nholes);

/**
* Create a geometry collection.
* \param type The geometry type, enumerated by \ref GEOSGeomTypes
* \param geoms A list of geometries that will form the collection
* \param ngeoms The number of geometries in the geoms list
* \return A newly allocated geometry collection. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \note The holes argument is an array of GEOSGeometry* objects.
*       The caller **retains ownership** of the containing array,
*       but the ownership of the pointed-to objects is transferred
*       to the returned \ref GEOSGeometry.
*/
GEOSGeometry* GEOSGeom_createCollection (
    int type,
    GEOSGeometry** geoms,
    uint ngeoms);

/**
* Create an empty geometry collection.
* \param type The geometry type, enumerated by \ref GEOSGeomTypes
* \return A newly allocated empty geometry collection. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_createEmptyCollection (int type);

/**
* Create a new copy of the input geometry.
* \param g The geometry to copy
* \return A newly allocated geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_clone (const(GEOSGeometry)* g);

/* ========== Memory management ========== */

/**
* Release the memory associated with a geometry.
* \param g The geometry to be destroyed.
*/
void GEOSGeom_destroy (GEOSGeometry* g);

/* ========== Topology Operations ========== */

/**
* Returns minimum rectangular polygon that contains the geometry.
* \param g The geometry to calculate an envelope for
* \return A newly allocated polygonal envelope. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSEnvelope (const(GEOSGeometry)* g);

/**
* Returns the intersection of two geometries: the set of points
* that fall within **both** geometries.
* \param g1 one of the geometries
* \param g2 the other geometry
* \return A newly allocated geometry of the intersection. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSIntersection (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* Returns the intersection of two geometries: the set of points
* that fall within **both** geometries. All the vertices of the output
* geometry must fall on the grid defined by the gridSize, and the
* output will be a valid geometry.
* \param g1 one of the geometries
* \param g2 the other geometry
* \param gridSize the cell size of the precision grid
* \return A newly allocated geometry of the intersection. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSIntersectionPrec (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2, double gridSize);

/**
* Returns the difference of two geometries A and B: the set of points
* that fall within A but **not** within B.
* \param ga the base geometry
* \param gb the geometry to subtract from it
* \return A newly allocated geometry of the difference. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSDifference (const(GEOSGeometry)* ga, const(GEOSGeometry)* gb);

/**
* Returns the difference of two geometries A and B: the set of points
* that fall within A but **not** within B.
* All the vertices of the output
* geometry must fall on the grid defined by the gridSize, and the
* output will be a valid geometry.
* \param ga one of the geometries
* \param gb the other geometry
* \param gridSize the cell size of the precision grid
* \return A newly allocated geometry of the difference. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSDifferencePrec (
    const(GEOSGeometry)* ga,
    const(GEOSGeometry)* gb,
    double gridSize);

/**
* Returns the symmetric difference of two geometries A and B: the set of points
* that fall in A but **not** within B and the set of points that fall in B but
* **not** in A.
* \param ga geometry A
* \param gb geometry B
* \return A newly allocated geometry of the symmetric difference. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSSymDifference (
    const(GEOSGeometry)* ga,
    const(GEOSGeometry)* gb);

/**
* Returns the symmetric difference of two geometries A and B: the set of points
* that fall in A but **not** within B and the set of points that fall in B but
* **not** in A.
* All the vertices of the output
* geometry must fall on the grid defined by the gridSize, and the
* output will be a valid geometry.
* \param ga one of the geometries
* \param gb the other geometry
* \param gridSize the cell size of the precision grid
* \return A newly allocated geometry of the symmetric difference. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSSymDifferencePrec (
    const(GEOSGeometry)* ga,
    const(GEOSGeometry)* gb,
    double gridSize);

/**
* Returns the union of two geometries A and B: the set of points
* that fall in A **or** within B.
* \param ga geometry A
* \param gb geometry B
* \return A newly allocated geometry of the union. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSUnion (const(GEOSGeometry)* ga, const(GEOSGeometry)* gb);

/**
* Returns the union of two geometries A and B: the set of points
* that fall in A **or** within B.
* All the vertices of the output
* geometry must fall on the grid defined by the gridSize, and the
* output will be a valid geometry.
* \param ga one of the geometries
* \param gb the other geometry
* \param gridSize the cell size of the precision grid
* \return A newly allocated geometry of the union. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSUnionPrec (
    const(GEOSGeometry)* ga,
    const(GEOSGeometry)* gb,
    double gridSize);

/**
* Returns the union of all components of a single geometry. Usually
* used to convert a collection into the smallest set of polygons
* that cover the same area.
* \param g The input geometry
* \return A newly allocated geometry of the union. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSUnaryUnion (const(GEOSGeometry)* g);

/**
* Returns the union of all components of a single geometry. Usually
* used to convert a collection into the smallest set of polygons
* that cover the same area.
* All the vertices of the output
* geometry must fall on the grid defined by the gridSize, and the
* output will be a valid geometry.
* \param g input geometry
* \param gridSize the cell size of the precision grid
* \return A newly allocated geometry of the union. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSUnaryUnionPrec (const(GEOSGeometry)* g, double gridSize);

/**
* Returns the "boundary" of a geometry, as defined by the DE9IM:
*
* - the boundary of a polygon is the linear rings dividing the exterior
*   from the interior
* - the boundary of a linestring is the end points
* - the boundary of a point is the point
*
* \param g The input geometry
* \return A newly allocated geometry of the boundary. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSBoundary (const(GEOSGeometry)* g);

/**
* Returns convex hull of a geometry. The smallest convex Geometry
* that contains all the points in the input Geometry
* \param g The input geometry
* \return A newly allocated geometry of the convex hull. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::overlayng::OverlayNG
*/
GEOSGeometry* GEOSConvexHull (const(GEOSGeometry)* g);

/**
* Returns the minimum rotated rectangular POLYGON which encloses
* the input geometry. The rectangle has width equal to the
* minimum diameter, and a longer length. If the convex hill of
* the input is degenerate (a line or point) a linestring or point
* is returned. The minimum rotated rectangle can be used as an
* extremely generalized representation for the given geometry.
* \param g The input geometry
* \return A newly allocated geometry of the rotated envelope. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSMinimumRotatedRectangle (const(GEOSGeometry)* g);

/**
* Constructs the "maximum inscribed circle" (MIC) for a polygonal geometry,
* up to a specified tolerance.
* The MIC is determined by a point in the interior of the area
* which has the farthest distance from the area boundary, along with a boundary point at that distance.
* In the context of geography the center of the MIC is known as the
* "pole of inaccessibility". A cartographic use case is to determine a suitable point
* to place a map label within a polygon.
* The radius length of the MIC is a  measure of how "narrow" a polygon is. It is the
* distance at which the negative buffer becomes empty.
* The class supports polygons with holes and multipolygons.
* The implementation uses a successive-approximation technique over a grid of square cells covering the area geometry.
* The grid is refined using a branch-and-bound algorithm. Point containment and distance are computed in a performant
* way by using spatial indexes.
* Returns a two-point linestring, with one point at the center of the inscribed circle and the other
* on the boundary of the inscribed circle.
* \param g Input geometry
* \param tolerance Stop the algorithm when the search area is smaller than this tolerance
* \return A newly allocated geometry of the MIC. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::algorithm::construct::MaximumInscribedCircle
*/
GEOSGeometry* GEOSMaximumInscribedCircle (
    const(GEOSGeometry)* g,
    double tolerance);

/**
* Constructs the "largest empty circle" (LEC) for a set of obstacle geometries, up to a
* specified tolerance. The obstacles are point and line geometries.
* The LEC is the largest circle which  has its center in the convex hull of the
* obstacles (the boundary), and whose interior does not intersect with any obstacle.
* The LEC center is the point in the interior of the boundary which has the farthest distance from
* the obstacles (up to tolerance). The LEC is determined by the center point and a point lying on an
* obstacle indicating the circle radius.
* The implementation uses a successive-approximation technique over a grid of square cells covering the obstacles and boundary.
* The grid is refined using a branch-and-bound algorithm.  Point containment and distance are computed in a performant
* way by using spatial indexes.
* Returns a two-point linestring, with one point at the center of the inscribed circle and the other
* on the boundary of the inscribed circle.
* \param obstacles The geometries that the LEC must fit within without covering
* \param boundary The area within which the LEC must reside
* \param tolerance Stop the algorithm when the search area is smaller than this tolerance
* \return A newly allocated geometry of the LEC. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::algorithm::construct::LargestEmptyCircle
*/
GEOSGeometry* GEOSLargestEmptyCircle (
    const(GEOSGeometry)* obstacles,
    const(GEOSGeometry)* boundary,
    double tolerance);

/**
* Returns a linestring geometry which represents the minimum diameter of the geometry.
* The minimum diameter is defined to be the width of the smallest band that
* contains the geometry, where a band is a strip of the plane defined
* by two parallel lines. This can be thought of as the smallest hole that the geometry
* can be moved through, with a single rotation.
* \param g The input geometry
* \return A newly allocated geometry of the LEC. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::algorithm::MinimumDiameter
*/
GEOSGeometry* GEOSMinimumWidth (const(GEOSGeometry)* g);

/**
* Computes the minimum clearance of a geometry.  The minimum clearance is the smallest amount by which
* a vertex could be move to produce an invalid polygon, a non-simple linestring, or a multipoint with
* repeated points.  If a geometry has a minimum clearance of 'eps', it can be said that:
*
* -  No two distinct vertices in the geometry are separated by less than 'eps'
* -  No vertex is closer than 'eps' to a line segment of which it is not an endpoint.
*
* If the minimum clearance cannot be defined for a geometry (such as with a single point, or a multipoint
* whose points are identical, a value of Infinity will be calculated.
*
* \param g the input geometry
* \param d a double to which the result can be stored
* \return 0 if no exception occurred.
*         2 if an exception occurred.
* \see geos::precision::MinimumClearance
*/
int GEOSMinimumClearance (const(GEOSGeometry)* g, double* d);

/**
* Returns a LineString whose endpoints define the minimum clearance of a geometry.
* If the geometry has no minimum clearance, an empty LineString will be returned.
*
* \param g the input geometry
* \return a linestring geometry, or NULL if an exception occurred.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::precision::MinimumClearance
*/
GEOSGeometry* GEOSMinimumClearanceLine (const(GEOSGeometry)* g);

/**
* Optimized union algorithm for polygonal inputs that are correctly
* noded and do not overlap. It will generate an error (return NULL)
* for inputs that do not satisfy this constraint.
* \param g The input geometry
* \return A geometry that covers all the points of the input geometry.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSCoverageUnion (const(GEOSGeometry)* g);

/**
* Returns a point that is inside the boundary of a polygonal geometry.
* \param g The input geometry
* \return A point that is inside the input
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::algorithm::InteriorPointArea
*/
GEOSGeometry* GEOSPointOnSurface (const(GEOSGeometry)* g);

/**
* Returns a point at the center of mass of the input.
* \param g The input geometry
* \return A point at the center of mass of the input
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::algorithm::Centroid
*/
GEOSGeometry* GEOSGetCentroid (const(GEOSGeometry)* g);

/**
* Returns a geometry which represents the "minimum bounding circle",
* the smallest circle that contains the input.
* \param[in] g The input geometry
* \param[out] radius Pointer will be filled with output radius.
* \param[out] center Pointer will be filled with output circle center. Caller must free.
* \return The circle geometry or NULL on exception
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::algorithm::MinimumBoundingCircle::getCircle
*/
GEOSGeometry* GEOSMinimumBoundingCircle (
    const(GEOSGeometry)* g,
    double* radius,
    GEOSGeometry** center);

/**
* For linear inputs, returns a new geometry in which no lines cross
* each other, and all touching occurs at end points.
* \param g The input geometry
* \return The noded geometry or NULL on exception
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::noding::GeometryNoder::node
*/
GEOSGeometry* GEOSNode (const(GEOSGeometry)* g);

/**
* Intersection optimized for a rectangular clipping polygon.
* Supposed to be faster than using GEOSIntersection(). Not
* guaranteed to return valid results.
* \param g The input geometry to be clipped
* \param xmin Left bound of clipping rectangle
* \param ymin Lower bound of clipping rectangle
* \param xmax Right bound of clipping rectangle
* \param ymax Upper bound of clipping rectangle
* \return The clipped geometry or NULL on exception
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::intersection::RectangleIntersection
*/
GEOSGeometry* GEOSClipByRect (
    const(GEOSGeometry)* g,
    double xmin,
    double ymin,
    double xmax,
    double ymax);

/**
* Polygonizes a set of Geometries which contain linework that
* represents the edges of a planar graph.
*
* All types of Geometry are accepted as input; the constituent
* linework is extracted as the edges to be polygonized.
*
* The edges must be correctly noded; that is, they must only meet
* at their endpoints. Polygonization will accept incorrectly noded
* input but will not form polygons from non-noded edges, and reports
* them as errors.
*
* The Polygonizer reports the following kinds of errors:
*
* - Dangles - edges which have one or both ends which are
*   not incident on another edge endpoint
* - Cut Edges - edges which are connected at both ends but
*   which do not form part of a polygon
* - Invalid Ring Lines - edges which form rings which are invalid
*   (e.g. the component lines contain a self-intersection)
*
* Errors are reported to output parameters "cuts", "dangles" and
* "invalid" (if not-null). Formed polygons are returned as a
* collection. NULL is returned on exception. All returned
* geometries must be destroyed by caller.
*
* The GEOSPolygonize_valid() variant allows extracting only polygons
* which form a valid polygonal result. The set of extracted polygons
* is guaranteed to be edge-disjoint. This is useful when it is known
* that the input lines form a valid polygonal geometry (which may
* include holes or nested polygons).
*
* \param geoms Array of linear geometries to polygons. Caller retains ownersihp of both array container and objects.
* \param ngeoms Size of the geoms array.
* \return The polygonal output geometry.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::polygonize::Polygonizer
*/
GEOSGeometry* GEOSPolygonize (const(GEOSGeometry*)* geoms, uint ngeoms);

/**
* Same polygonizing behavior as GEOSPolygonize(), but only returning results
* that are valid.
*
* \param geoms Array of linear geometries to polygons. Caller retains ownersihp of both array container and objects.
* \param ngeoms Size of the geoms array.
* \return The polygonal output geometry.
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::polygonize::Polygonizer
*/
GEOSGeometry* GEOSPolygonize_valid (const(GEOSGeometry*)* geoms, uint ngeoms);

/**
* Perform the polygonization as GEOSPolygonize() but return only the
* "cut edges", the linear features that are connected at both ends,
* do *not* participate in the final polygon.
*
* \param geoms Array of linear geometries to polygons. Caller retains ownersihp of both array container and objects.
* \param ngeoms Size of the geoms array.
* \return The "cut edges"
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::polygonize::Polygonizer
*/
GEOSGeometry* GEOSPolygonizer_getCutEdges (
    const(GEOSGeometry*)* geoms,
    uint ngeoms);

/**
* Perform the polygonization as GEOSPolygonize() and return the
* polygonal result as well as all extra ouputs.
*
* \param[in] input A single geometry with all the input lines to polygonize.
* \param[out] cuts Pointer to hold "cut edges", connected on both ends but not part of output. Caller must free.
* \param[out] dangles Pointer to hold "dangles", connected one end but not part of output. Caller must free.
* \param[out] invalid Pointer to hold invalid outputs, polygons formed but not valid. Caller must free.
* \return The polygonal valid output
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::polygonize::Polygonizer
*/
GEOSGeometry* GEOSPolygonize_full (
    const(GEOSGeometry)* input,
    GEOSGeometry** cuts,
    GEOSGeometry** dangles,
    GEOSGeometry** invalid);

/**
* Perform a polygonization using all the linework, assuming that
* rings contained within rings are empty holes, rather then
* extra polygons.
* \param g The input linework
* \return The polygonal output
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::polygonize::BuildArea
*/
GEOSGeometry* GEOSBuildArea (const(GEOSGeometry)* g);

/**
* Sews together a set of fully noded LineStrings
* removing any cardinality 2 nodes in the linework.
* \param g The input linework
* \return The merged linework
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::linemerge::LineMerger
*/
GEOSGeometry* GEOSLineMerge (const(GEOSGeometry)* g);

/**
* For geometries with coordinate sequences, reverses the order
* of the sequences. Converts CCW rings to CW. Reverses direction
* of LineStrings.
* \param g The input geometry
* \return The reversed geometry
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSReverse (const(GEOSGeometry)* g);

/**
* Apply the
* [Douglas/Peucker algorithm](https://en.wikipedia.org/wiki/Ramer–Douglas–Peucker_algorithm)
* to the coordinate sequences of the input geometry.
* Removes "unnecessary" vertices, vertices
* that are co-linear within the tolerance distance.
* \param g The input geometry
* \param tolerance The tolerance to apply. Larger tolerance leads to simpler output.
* \return The simplified geometry
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::simplify::DouglasPeuckerSimplifier
*/
GEOSGeometry* GEOSSimplify (const(GEOSGeometry)* g, double tolerance);

/**
* Apply the
* [Douglas/Peucker algorithm](https://en.wikipedia.org/wiki/Ramer–Douglas–Peucker_algorithm)
* to the coordinate sequences of the input geometry.
* Removes "unnecessary" vertices, vertices
* that are co-linear within the tolerance distance.
* Returns a valid output geometry, checking for collapses, ring-intersections, etc
* and attempting to avoid. More computationally expensive than GEOSSimplify()
* \param g The input geometry
* \param tolerance The tolerance to apply. Larger tolerance leads to simpler output.
* \return The simplified geometry
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::simplify::DouglasPeuckerSimplifier
*/
GEOSGeometry* GEOSTopologyPreserveSimplify (
    const(GEOSGeometry)* g,
    double tolerance);

/**
* Return all distinct vertices of input geometry as a MultiPoint.
* Note that only 2 dimensions of the vertices are considered when
* testing for equality.
* \param g The input geometry
* \return The distinct points
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSGeom_extractUniquePoints (const(GEOSGeometry)* g);

/**
* Find paths shared between the two given lineal geometries.
*
* Returns a GeometryCollection having two elements:
*
* - first element is a MultiLineString containing shared paths
*   having the _same_ direction on both inputs
* - second element is a MultiLineString containing shared paths
*   having the _opposite_ direction on the two inputs
*
* \param g1 An input geometry
* \param g2 An input geometry
* \return The shared paths
* Caller is responsible for freeing with GEOSGeom_destroy().
* \see geos::operation::sharedpaths::SharedPathsOp
*/
GEOSGeometry* GEOSSharedPaths (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/**
* Snap first geometry onto second within the given tolerance.
* \param input An input geometry
* \param snap_target A geometry to snap the input to
* \param tolerance Snapping tolerance
* \return The snapped verion of the input. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSSnap (
    const(GEOSGeometry)* input,
    const(GEOSGeometry)* snap_target,
    double tolerance);

/**
* Return a Delaunay triangulation of the vertices of the given geometry.
*
* \param g the input geometry whose vertices will be used as "sites"
* \param tolerance optional snapping tolerance to use for improved robustness
* \param onlyEdges if non-zero will return a MultiLineString, otherwise it will
*                  return a GeometryCollection containing triangular Polygons.
*
* \return A newly allocated geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSDelaunayTriangulation (
    const(GEOSGeometry)* g,
    double tolerance,
    int onlyEdges);

/**
* Return a constrained Delaunay triangulation of the vertices of the
* given polygon(s).
* For non-polygonal inputs, returns an empty geometry collection.
*
* \param g the input geometry whose rings will be used as input
* \return A newly allocated geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSConstrainedDelaunayTriangulation (const(GEOSGeometry)* g);

/**
* Returns the Voronoi polygons of the vertices of the given geometry.
*
* \param g the input geometry whose vertices will be used as sites.
* \param tolerance snapping tolerance to use for improved robustness
* \param onlyEdges whether to return only edges of the voronoi cells
* \param env clipping envelope for the returned diagram, automatically
*            determined if env is NULL.
*            The diagram will be clipped to the larger
*            of this envelope or an envelope surrounding the sites.
*
* \return A newly allocated geometry. NULL on exception.
* Caller is responsible for freeing with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSVoronoiDiagram (
    const(GEOSGeometry)* g,
    const(GEOSGeometry)* env,
    double tolerance,
    int onlyEdges);

/**
* Computes the coordinate where two line segments intersect, if any
*
* \param[in] ax0 x-coordinate of 1st point in 1st segment
* \param[in] ay0 y-coordinate of 1st point in 1st segment
* \param[in] ax1 x-coordinate of 2nd point in 1st segment
* \param[in] ay1 y-coordinate of 2nd point in 1st segment
* \param[in] bx0 x-coordinate of 1st point in 2nd segment
* \param[in] by0 y-coordinate of 1st point in 2nd segment
* \param[in] bx1 x-coordinate of 2nd point in 2nd segment
* \param[in] by1 y-coordinate of 2nd point in 2nd segment
* \param[out] cx x-coordinate of intersection point
* \param[out] cy y-coordinate of intersection point
*
* \return 0 on error, 1 on success, -1 if segments do not intersect
*/
int GEOSSegmentIntersection (
    double ax0,
    double ay0,
    double ax1,
    double ay1,
    double bx0,
    double by0,
    double bx1,
    double by1,
    double* cx,
    double* cy);

/* ========= Binary Predicates ========= */

/**
* True if no point of either geometry touchess or is within the other.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::disjoint
*/
char GEOSDisjoint (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometries share boundaries at one or more points, but do
* not have interior overlaps.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::touches
*/
char GEOSTouches (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometries are not disjoint.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::intersects
*/
char GEOSIntersects (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometries interiors interact but their boundares do not.
* Most useful for finding line crosses cases.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::crosses
*/
char GEOSCrosses (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometry g1 is completely within g2, and not
* touching the boundary of g2.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::within
*/
char GEOSWithin (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometry g2 is completely within g1.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::contains
*/
char GEOSContains (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometries share interiors but are neither
* within nor contained.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::overlaps
*/
char GEOSOverlaps (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometries cover the same space on the place.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::equals
*/
char GEOSEquals (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometry g1 is completely within g2, including possibly
* touching the boundary of g2.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::covers
*/
char GEOSCovers (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* True if geometry g2 is completely within g1, including possibly
* touching the boundary of g1.
* \param g1 Input geometry
* \param g2 Input geometry
* \returns 1 on true, 0 on false, 2 on exception
* \see geos::geom::Geometry::coveredby
*/
char GEOSCoveredBy (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* Determine pointwise equivalence of two geometries, by checking if each vertex of g2 is
* within tolerance of the corresponding vertex in g1.
* Unlike GEOSEquals(), geometries that are topologically equivalent but have different
* representations (e.g., LINESTRING (0 0, 1 1) and MULTILINESTRING ((0 0, 1 1)) ) are not
* considered equivalent by GEOSEqualsExact().
* \param g1 Input geometry
* \param g2 Input geometry
* \param tolerance Tolerance to determine vertex equality
* \returns 1 on true, 0 on false, 2 on exception
*/
char GEOSEqualsExact (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double tolerance);

/* ========== Prepared Geometry Binary predicates ========== */

/**
* A \ref GEOSPreparedGeometry is a wrapper around \ref GEOSGeometry
* that adds in a spatial index on the edges of the geometry. This
* internal index allows spatial predicates to evaluate much faster,
* so for cases in which the same base geometry will be used over and
* over again for predicate tests, wrapping it in a \ref GEOSPreparedGeometry
* is a best practice.
*
* The caller retains ownership of the base geometry, and after
* processing is complete, must destroy **both** the prepared and the
* base geometry. (Ideally, destroy the prepared geometry first, as
* it has an internal reference to the base geometry.)
*
* \param g The base geometry to wrap in a prepared geometry.
* \return A prepared geometry. Caller is responsible for freeing with
*         GEOSPreparedGeom_destroy()
*/
const(GEOSPreparedGeometry)* GEOSPrepare (const(GEOSGeometry)* g);

/**
* Free the memory associated with a \ref GEOSPreparedGeometry.
* Caller must separately free the base \ref GEOSGeometry used
* to create the prepared geometry.
* \param g Prepared geometry to destroy.
*/
void GEOSPreparedGeom_destroy (const(GEOSPreparedGeometry)* g);

/**
* Using a \ref GEOSPreparedGeometry do a high performance
* calculation of whether the provided geometry is contained.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSContains
*/
char GEOSPreparedContains (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedGeometry do a high performance
* calculation of whether the provided geometry is contained properly.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSContainsProperly
*/
char GEOSPreparedContainsProperly (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedGeometry do a high performance
* calculation of whether the provided geometry is covered by.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSCoveredBy
*/
char GEOSPreparedCoveredBy (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedGeometry do a high performance
* calculation of whether the provided geometry covers.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSCovers
*/
char GEOSPreparedCovers (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedGeometry do a high performance
* calculation of whether the provided geometry crosses.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSCrosses
*/
char GEOSPreparedCrosses (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedDisjoint do a high performance
* calculation of whether the provided geometry is disjoint.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSDisjoin
*/
char GEOSPreparedDisjoint (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedDisjoint do a high performance
* calculation of whether the provided geometry is disjoint.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSDisjoin
*/
char GEOSPreparedIntersects (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedDisjoint do a high performance
* calculation of whether the provided geometry overlaps.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSOverlaps
*/
char GEOSPreparedOverlaps (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedDisjoint do a high performance
* calculation of whether the provided geometry touches.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSTouches
*/
char GEOSPreparedTouches (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedDisjoint do a high performance
* calculation of whether the provided geometry is within.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns 1 on true, 0 on false, 2 on exception
* \see GEOSWithin
*/
char GEOSPreparedWithin (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedDisjoint do a high performance
* calculation to find the nearest points between the
* prepared and provided geometry.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \returns A coordinate sequence containing the nearest points, or NULL on exception.
*          The first point in the sequence is from the prepared geometry, and the
*          seconds is from the other argument.
*/
GEOSCoordSequence* GEOSPreparedNearestPoints (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2);

/**
* Using a \ref GEOSPreparedDistance do a high performance
* calculation to find the distance points between the
* prepared and provided geometry. Useful for situations where
* one geometry is large and static and needs to be tested
* against a large number of other geometries.
* \param[in] pg1 The prepared geometry
* \param[in] g2 The geometry to test
* \param[out] dist Pointer to store the result in
* \return 1 on success
*/
int GEOSPreparedDistance (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2,
    double* dist);

/**
* Using a \ref GEOSPreparedDistanceWithin do a high performance
* calculation to find whether the prepared and provided geometry
* are within the given max distance.
* Useful for situations where
* one geometry is large and static and needs to be tested
* against a large number of other geometries.
* \param pg1 The prepared geometry
* \param g2 The geometry to test
* \param dist The max distance
* \return 1 on success
*/
char GEOSPreparedDistanceWithin (
    const(GEOSPreparedGeometry)* pg1,
    const(GEOSGeometry)* g2,
    double dist);

/* ========== STRtree functions ========== */

/**
* Create a new \ref GEOSSTRtree using the Sort-Tile-Recursive algorithm
* ([STRtree](https://en.wikipedia.org/wiki/R-tree))
* for two-dimensional spatial data.
*
* \param nodeCapacity The maximum number of child nodes that a node may have.
*        The minimum recommended capacity value is 4.
*        If unsure, use a default node capacity of 10.
* \return a pointer to the created tree
*/
GEOSSTRtree* GEOSSTRtree_create (size_t nodeCapacity);

/**
* Insert an item into an \ref GEOSSTRtree
*
* \param tree the \ref GEOSSTRtree in which the item should be inserted
* \param g a GEOSGeometry whose envelope corresponds to the extent of 'item'. As of GEOS 3.9, this envelope will be
 *       copied into the tree and the caller may destroy `g` while the tree is still in use. Before GEOS 3.9, `g`
 *       must be retained until the tree is destroyed.
* \param item the item to insert into the tree
* \note The tree does **not** take ownership of the geometry or the item.
*/
void GEOSSTRtree_insert (GEOSSTRtree* tree, const(GEOSGeometry)* g, void* item);

/**
* Query an \ref GEOSSTRtree for items intersecting a specified envelope
*
* \param tree the \ref GEOSSTRtree to search
* \param g a GEOSGeomety from which a query envelope will be extracted
* \param callback a function to be executed for each item in the tree whose envelope intersects
*            the envelope of 'g'.  The callback function should take two parameters: a void
*            pointer representing the located item in the tree, and a void userdata pointer.
* \param userdata an optional pointer to pe passed to 'callback' as an argument
*/
void GEOSSTRtree_query (
    GEOSSTRtree* tree,
    const(GEOSGeometry)* g,
    GEOSQueryCallback callback,
    void* userdata);

/**
* Returns the nearest item in the \ref GEOSSTRtree to the supplied geometry.
* All items in the tree MUST be of type \ref GEOSGeometry.
* If this is not the case, use GEOSSTRtree_nearest_generic() instead.
*
* \param tree the \ref GEOSSTRtree to search
* \param geom the geometry with which the tree should be queried
* \return a const pointer to the nearest \ref GEOSGeometry in the tree to 'geom', or NULL in
*            case of exception
*/
const(GEOSGeometry)* GEOSSTRtree_nearest (
    GEOSSTRtree* tree,
    const(GEOSGeometry)* geom);

/**
* Returns the nearest item in the \ref GEOSSTRtree to the supplied item
*
* \param tree the STRtree to search
* \param item the item with which the tree should be queried
* \param itemEnvelope a GEOSGeometry having the bounding box of 'item'
* \param distancefn a function that can compute the distance between two items
*            in the STRtree.  The function should return zero in case of error,
*            and should store the computed distance to the location pointed to by
*            the 'distance' argument.  The computed distance between two items
*            must not exceed the Cartesian distance between their envelopes.
* \param userdata optional pointer to arbitrary data; will be passed to distancefn
*            each time it is called.
* \return a const pointer to the nearest item in the tree to 'item', or NULL in
*            case of exception
*/
const(void)* GEOSSTRtree_nearest_generic (
    GEOSSTRtree* tree,
    const(void)* item,
    const(GEOSGeometry)* itemEnvelope,
    GEOSDistanceCallback distancefn,
    void* userdata);

/**
* Iterate over all items in the \ref GEOSSTRtree.
*
* \param tree the STRtree over which to iterate
* \param callback a function to be executed for each item in the tree.
* \param userdata payload to pass the callback function.
*/
void GEOSSTRtree_iterate (
    GEOSSTRtree* tree,
    GEOSQueryCallback callback,
    void* userdata);

/**
 * Removes an item from the \ref GEOSSTRtree
 *
 * \param tree the STRtree from which to remove an item
 * \param g the envelope of the item to remove
 * \param item the item to remove
 * \return 0 if the item was not removed;
 *         1 if the item was removed;
 *         2 if an exception occurred
 */
char GEOSSTRtree_remove (GEOSSTRtree* tree, const(GEOSGeometry)* g, void* item);

/**
* Frees all the memory associated with a \ref GEOSSTRtree.
* Only the tree is freed. The geometries and items fed into
* GEOSSTRtree_insert() are not owned by the tree, and are
* still left to the caller to manage.
*/
void GEOSSTRtree_destroy (GEOSSTRtree* tree);

/* ========= Unary predicates ========= */

/**
* Tests whether the input geometry is empty. If the geometry or any
* component is non-empty, the geometry is non-empty. An empty geometry
* has no boundary or interior.
* \param g The geometry to test
* \return 1 on true, 0 on false, 2 on exception
*/
char GEOSisEmpty (const(GEOSGeometry)* g);

/**
* Tests whether the input geometry is "simple". Mostly relevant for
* linestrings. A "simple" linestring has no self-intersections.
* \param g The geometry to test
* \return 1 on true, 0 on false, 2 on exception
*/
char GEOSisSimple (const(GEOSGeometry)* g);

/**
* Tests whether the input geometry is a ring. Rings are
* linestrings, without self-intersections,
* with start and end point being identical.
* \param g The geometry to test
* \return 1 on true, 0 on false, 2 on exception
*/
char GEOSisRing (const(GEOSGeometry)* g);

/**
* Tests whether the input geometry has z coordinates.
* \param g The geometry to test
* \return 1 on true, 0 on false, 2 on exception
*/
char GEOSHasZ (const(GEOSGeometry)* g);

/**
* Tests whether the input geometry is closed.
* A closed geometry is a linestring or multilinestring
* with the start and end points being the same.
* \param g The geometry to test
* \return 1 on true, 0 on false, 2 on exception
*/
char GEOSisClosed (const(GEOSGeometry)* g);

/* ========= Dimensionally Extended 9 Intersection Model ========= */

/**
* Calculate the DE9IM pattern for this geometry pair
* and compare against the provided pattern to check for
* consistency. If the result and pattern are consistent
* return true. The pattern may include glob "*" characters
* for portions that are allowed to match any value.
* \see geos::geom::Geometry::relate
* \param g1 First geometry in pair
* \param g2 Second geometry in pair
* \param pat DE9IM pattern to check
* \return 1 on true, 0 on false, 2 on exception
*/
char GEOSRelatePattern (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    const(char)* pat);

/**
* Calculate and return the DE9IM pattern for this geometry pair.
* \see geos::geom::Geometry::relate
* \param g1 First geometry in pair
* \param g2 Second geometry in pair
* \return DE9IM string. Caller is responsible for freeing with GEOSFree().
*         NULL on exception
*/
char* GEOSRelate (const(GEOSGeometry)* g1, const(GEOSGeometry)* g2);

/**
* Compare two DE9IM patterns and return true if they
* are consistent.
* \param mat Complete DE9IM string (does not have "*")
* \param pat Pattern to match to (may contain "*")
* \return 1 on true, 0 on false, 2 on exception
*/
char GEOSRelatePatternMatch (const(char)* mat, const(char)* pat);

/**
* Calculate and return the DE9IM pattern for this geometry pair.
* Apply the supplied \ref GEOSRelateBoundaryNodeRules.
* \see geos::geom::Geometry::relate
* \see geos::algorithm::BoundaryNodeRule
* \param g1 First geometry in pair
* \param g2 Second geometry in pair
* \param bnr A member of the \ref GEOSRelateBoundaryNodeRules enum
* \return DE9IM string. Caller is responsible for freeing with GEOSFree().
*         NULL on exception
*/
char* GEOSRelateBoundaryNodeRule (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    int bnr);

/* ========== Validity checking ========== */

/**
* Check the validity of the provided geometry.
* - All points are valid.
* - All non-zero-length linestrings are valid.
* - Polygon rings must be non-self-intersecting, and interior rings
*   contained within exterior rings.
* - Multi-polygon components may not touch or overlap.
*
* \param g The geometry to test
* \return 1 on true, 0 on false, 2 on exception
* \see geos::operation::valid::isValidOp
*/
char GEOSisValid (const(GEOSGeometry)* g);

/**
* Return the human readable reason a geometry is invalid,
* "Valid Geometry" string otherwise, or NULL on exception.
* \param g The geometry to test
* \return A string with the reason, NULL on exception.
          Caller must GEOSFree() their result.
*/
char* GEOSisValidReason (const(GEOSGeometry)* g);

/**
* In one step, calculate and return the validity, the
* human readable validity reason and a point at which validity
* rules are broken.
* Caller has the responsibility to destroy 'reason' with GEOSFree()
* and 'location' with GEOSGeom_destroy()
* \param g The geometry to test
* \param flags A value from the \ref GEOSValidFlags enum
* \param reason A pointer in which the reason string will be places
* \param location A pointer in which the location GEOSGeometry will be placed
* \return 1 when valid, 0 when invalid, 2 on exception
*/
char GEOSisValidDetail (
    const(GEOSGeometry)* g,
    int flags,
    char** reason,
    GEOSGeometry** location);

/**
* Repair an invalid geometry, returning a valid output.
* \param g The geometry to repair
* \return The repaired geometry. Caller must free with GEOSGeom_destroy().
*/
GEOSGeometry* GEOSMakeValid (const(GEOSGeometry)* g);

/**
* Repair an invalid geometry, returning a valid output, using the
* indicated GEOSMakeValidMethods algorithm and options.
* \param g is the geometry to test.
* \param makeValidParams is a GEOSMakeValidParams with the desired parameters set on it.
* \return A repaired geometry. Caller must free with GEOSGeom_destroy().
* \see GEOSMakeValidParams_create
* \see GEOSMakeValidParams_destroy
* \see GEOSMakeValidParams_setMethod
* \see GEOSMakeValidParams_setKeepCollapsed
*/
GEOSGeometry* GEOSMakeValidWithParams (
    const(GEOSGeometry)* g,
    const(GEOSMakeValidParams)* makeValidParams);

/**
* Create a GEOSMakeValidParams to hold the desired parameters
* to control the algorithm and behavior of the validation process.
* \return a parameter object
* \see GEOSMakeValidWithParams
*/
GEOSMakeValidParams* GEOSMakeValidParams_create ();

/**
* Destroy a GEOSMakeValidParams.
* \param parms the object to destroy
* \see GEOSMakeValidWithParams
*/
void GEOSMakeValidParams_destroy (GEOSMakeValidParams* parms);

/**
* Set the GEOSMakeValidMethods to use in making the geometry
* valid.
* \return 0 on exception, 1 on success.
* \see GEOSMakeValidWithParams
*/
int GEOSMakeValidParams_setMethod (
    GEOSMakeValidParams* p,
    GEOSMakeValidMethods method);

/**
* When this parameter is not set to 0, the GEOS_MAKE_VALID_STRUCTURE method will drop
* any component that has collapsed into a lower dimensionality.
* For example, a ring collapsing to a line, or a line collapsing
* to a point.
* \return 0 on exception, 1 on success.
* \see GEOSMakeValidWithParams
*/
int GEOSMakeValidParams_setKeepCollapsed (
    GEOSMakeValidParams* p,
    int keepCollapsed);

/* ========== Geometry info ========== */

/**
* Returns the geometry type string for this geometry.
* eg: "GeometryCollection", "LineString"
* \param g Input geometry
* \return A string with the geometry type.
* Caller must free with GEOSFree().
* NULL on exception.
*/
char* GEOSGeomType (const(GEOSGeometry)* g);

/**
* Returns the \ref GEOSGeomTypeId number for this geometry.
* \param g Input geometry
* \return The geometry type number, or -1 on exception.
*/
int GEOSGeomTypeId (const(GEOSGeometry)* g);

/**
* Returns the "spatial reference id" (SRID) for this geometry.
* \param g Input geometry
* \return SRID number or 0 if unknown / not set.
*/
int GEOSGetSRID (const(GEOSGeometry)* g);

/**
* Set the "spatial reference id" (SRID) for this geometry.
* \param g Input geometry
* \param SRID SRID number or 0 for unknown SRID.
*/
void GEOSSetSRID (GEOSGeometry* g, int SRID);

/**
* Return the anonymous "user data" for this geometry.
* User data must be managed by the caller, and freed before
* the geometry is freed.
* \param g Input geometry
* \return A void* to the user data, caller is responsible for
*         casting to the appropriate type and freeing.
*/
void* GEOSGeom_getUserData (const(GEOSGeometry)* g);

/**
* Set the anonymous "user data" for this geometry.
* Don't forget to free the user data before freeing the geometry.
* \param g Input geometry
* \param userData Void pointer to user data
*/
void GEOSGeom_setUserData (GEOSGeometry* g, void* userData);

/**
* Returns the number of sub-geometries immediately under a
* multi-geometry or collection or 1 for a simple geometry.
* For nested collections, remember to check if returned
* sub-geometries are **themselves** also collections.
* \param g Input geometry
* \return Number of direct children in this collection
* \warning For GEOS < 3.2 this function may crash when fed simple geometries
*/
int GEOSGetNumGeometries (const(GEOSGeometry)* g);

/**
* Returns the specified sub-geometry of a collection. For
* a simple geometry, returns a pointer to the input.
* Returned object is a pointer to internal storage:
* it must NOT be destroyed directly.
* \param g Input geometry
* \param n Sub-geometry index, zero-base
* \return A const \ref GEOSGeometry, do not free!
          It will be freed when the parent is freed.
          Returns NULL on exception.
* \note Up to GEOS 3.2.0 the input geometry must be a Collection, in
*       later versions it doesn't matter (getGeometryN(0) for a single will
*       return the input).
*/
const(GEOSGeometry)* GEOSGetGeometryN (const(GEOSGeometry)* g, int n);

/**
* Organize the sub-geometries, rings, and coordinate order
* of geometries, so that geometries that represent the same
* object can be easily compared. Starts rings from the same
* location, orients them in the same way, sorts geometry
* sub-components in the same way. Use before calling
* GEOSEqualsExact() to avoid false "not equal" results.
* \param g Input geometry
* \return 0 on success or -1 on exception
*/
int GEOSNormalize (GEOSGeometry* g);

/**
* Change the rounding precision on a geometry. This will
* affect the precision of the existing geometry as well as
* any geometries derived from this geometry using overlay
* functions. The output will be a valid \ref GEOSGeometry.
*
* Note that operations will always be performed in the precision
* of the geometry with higher precision (smaller "gridSize").
* That same precision will be attached to the operation outputs.
*
* In the Default and GEOS_PREC_KEEP_COLLAPSED modes invalid input
* may cause an error to occur, unless the invalidity is below
* the scale of the requested precision
*
* There are only 3 modes. The GEOS_PREC_NO_TOPO mode
* takes precedence over GEOS_PREC_KEEP_COLLAPSED.
* So the combination GEOS_PREC_NO_TOPO || GEOS_PREC_KEEP_COLLAPSED
* has the same semantics as GEOS_PREC_NO_TOPO
*
* \param g Input geometry
* \param gridSize cell size of grid to round coordinates to,
*        or 0 for FLOATING precision
* \param flags The bitwise OR of members of the \ref GEOSPrecisionRules enum
* \return The precision reduced result.
*         Caller must free with GEOSGeom_destroy()
*         NULL on exception.
*/
GEOSGeometry* GEOSGeom_setPrecision (
    const(GEOSGeometry)* g,
    double gridSize,
    int flags);

/**
* Read the currently set precision value from the
* geometry and returns the grid size if it is a fixed
* precision or 0.0 if it is full floating point precision.
* \param g Input geometry
* \return The grid size, or -1 on exception
*/
double GEOSGeom_getPrecision (const(GEOSGeometry)* g);

/**
* Returns the number of interior rings, for a Polygon input, or
* an exception otherwise.
* \param g Input Polygon geometry
* \return Number of interior rings, -1 on exception
*/
int GEOSGetNumInteriorRings (const(GEOSGeometry)* g);

/**
* Returns the number of points, for a LineString input, or
* an exception otherwise.
* \param g Input LineString geometry
* \return Number of points, -1 on exception
*/
int GEOSGeomGetNumPoints (const(GEOSGeometry)* g);

/**
* Returns the X coordinate, for a Point input, or an
* exception otherwise.
* \param[in] g Input Point geometry
* \param[out] x Pointer to hold return value
* \returns 1 on success, 0 on exception
*/
int GEOSGeomGetX (const(GEOSGeometry)* g, double* x);

/**
* Returns the Y coordinate, for a Point input, or an
* exception otherwise.
* \param[in] g Input Point geometry
* \param[out] y Pointer to hold return value
* \returns 1 on success, 0 on exception
*/
int GEOSGeomGetY (const(GEOSGeometry)* g, double* y);

/**
* Returns the Z coordinate, for a Point input, or an
* exception otherwise.
* \param[in] g Input Point geometry
* \param[out] z Pointer to hold return value
* \returns 1 on success, 0 on exception
*/
int GEOSGeomGetZ (const(GEOSGeometry)* g, double* z);

/**
* Returns the N'th ring for a Polygon input.
* \note Returned object is a pointer to internal storage:
*       it must NOT be destroyed directly.
* \param g Input Polygon geometry
* \param n Index of the desired ring
* \return LinearRing geometry. Owned by parent geometry, do not free. NULL on exception.
*/
const(GEOSGeometry)* GEOSGetInteriorRingN (const(GEOSGeometry)* g, int n);

/**
* Get the external ring of a Polygon.
* \note Returned object is a pointer to internal storage:
*       it must NOT be destroyed directly.
* \param g Input Polygon geometry
* \return LinearRing geometry. Owned by parent geometry, do not free. NULL on exception.
*/
const(GEOSGeometry)* GEOSGetExteriorRing (const(GEOSGeometry)* g);

/**
* Get the total number of points in a geometry,
* of any type.
* \param g Input geometry
* \return Number of points in the geometry. -1 on exception.
*/
int GEOSGetNumCoordinates (const(GEOSGeometry)* g);

/**
* Return the coordinate sequence underlying the
* given geometry (Must be a LineString, LinearRing or Point).
* Do not directly free the coordinate sequence, it is owned by
* the parent geometry.
* \param g Input geometry
* \return Coordinate sequence or NULL on exception.
*/
const(GEOSCoordSequence)* GEOSGeom_getCoordSeq (const(GEOSGeometry)* g);

/**
* Return the planar dimensionality of the geometry.
*
* - 0 for point, multipoint
* - 1 for linestring, multilinestring
* - 2 for polygon, multipolygon
*
* \see geos::geom::Dimension::DimensionType
* \param g Input geometry
* \return The dimensionality
*/
int GEOSGeom_getDimensions (const(GEOSGeometry)* g);

/**
* Return the cartesian dimension of the geometry.
*
* - 2 for XY data
* - 3 for XYZ data
*
* \param g Input geometry
* \return The dimension
*/
int GEOSGeom_getCoordinateDimension (const(GEOSGeometry)* g);

/**
* Finds the minimum X value in the geometry.
* \param[in] g Input geometry
* \param[out] value Pointer to place result
* \return 0 on exception
*/
int GEOSGeom_getXMin (const(GEOSGeometry)* g, double* value);

/**
* Finds the minimum Y value in the geometry.
* \param[in] g Input geometry
* \param[out] value Pointer to place result
* \return 0 on exception
*/
int GEOSGeom_getYMin (const(GEOSGeometry)* g, double* value);

/**
* Finds the maximum X value in the geometry.
* \param[in] g Input geometry
* \param[out] value Pointer to place result
* \return 0 on exception
*/
int GEOSGeom_getXMax (const(GEOSGeometry)* g, double* value);

/**
* Finds the maximum Y value in the geometry.
* \param[in] g Input geometry
* \param[out] value Pointer to place result
* \return 0 on exception
*/
int GEOSGeom_getYMax (const(GEOSGeometry)* g, double* value);

/**
* Return the N'th point of a LineString
* \param g Input geometry, must be a LineString
* \param n Index of desired point (zero based)
* \return A Point geometry.
*         Caller must free with GEOSGeom_destroy()
*         NULL on exception.
*/
GEOSGeometry* GEOSGeomGetPointN (const(GEOSGeometry)* g, int n);

/**
* Return the first point of a LineString
* \param g Input geometry, must be a LineString
* \return A Point geometry.
*         Caller must free with GEOSGeom_destroy()
*         NULL on exception.
*/
GEOSGeometry* GEOSGeomGetStartPoint (const(GEOSGeometry)* g);

/**
* Return the last point of a LineString
* \param g Input geometry, must be a LineString
* \return A Point geometry.
*         Caller must free with GEOSGeom_destroy()
*         NULL on exception.
*/
GEOSGeometry* GEOSGeomGetEndPoint (const(GEOSGeometry)* g);

/* ========= Misc functions ========= */

/**
* Calculate the area of a geometry.
* \param[in] g Input geometry
* \param[out] area Pointer to be filled in with area result
* \return 1 on success, 0 on exception.
*/
int GEOSArea (const(GEOSGeometry)* g, double* area);

/**
* Calculate the length of a geometry.
* \param[in] g Input geometry
* \param[out] length Pointer to be filled in with length result
* \return 1 on success, 0 on exception.
*/
int GEOSLength (const(GEOSGeometry)* g, double* length);

/**
* Calculate the distance between two geometries.
* \param[in] g1 Input geometry
* \param[in] g2 Input geometry
* \param[out] dist Pointer to be filled in with distance result
* \return 1 on success, 0 on exception.
*/
int GEOSDistance (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/**
* Test whether the distance between two geometries is
* within the given dist.
* \param g1 Input geometry
* \param g2 Input geometry
* \param dist The max distance
* \returns 1 on true, 0 on false, 2 on exception
*/
char GEOSDistanceWithin (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double dist);

/**
* Calculate the distance between two geometries, using the
* indexed facet distance, which first indexes the geometries
* internally, then calculates the distance. Useful when one
* or both geometries is very large.
* \param[in] g1 Input geometry
* \param[in] g2 Input geometry
* \param[out] dist Pointer to be filled in with distance result
* \return 1 on success, 0 on exception.
* \see geos::operation::distance:;IndexedFacetDistance
*/
int GEOSDistanceIndexed (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/**
* Calculate the Hausdorff distance between two geometries.
* [Hausdorff distance](https://en.wikipedia.org/wiki/Hausdorff_distance)
* is the largest distance between two geometries.
* \param[in] g1 Input geometry
* \param[in] g2 Input geometry
* \param[out] dist Pointer to be filled in with distance result
* \return 1 on success, 0 on exception.
* \see geos::algorithm::distance::DiscreteHausdorffDistance
*/
int GEOSHausdorffDistance (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/**
* Calculate a more precise Hausdorff distance between two geometries,
* by densifying the inputs before computation.
* [Hausdorff distance](https://en.wikipedia.org/wiki/Hausdorff_distance)
* is the largest distance between two geometries.
* \param[in] g1 Input geometry
* \param[in] g2 Input geometry
* \param[in] densifyFrac The largest % of the overall line length that
*            any given two-point segment should be
* \param[out] dist Pointer to be filled in with distance result
* \return 1 on success, 0 on exception.
* \see geos::algorithm::distance::DiscreteHausdorffDistance
*/
int GEOSHausdorffDistanceDensify (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double densifyFrac,
    double* dist);

/**
* Calculate the
* [Frechet distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance)
* between two geometries,
* a similarity measure for linear features.
* \param[in] g1 Input geometry
* \param[in] g2 Input geometry
* \param[out] dist Pointer to be filled in with distance result
* \return 1 on success, 0 on exception.
* \see geos::algorithm::distance::DiscreteFrechetDistance
*/
int GEOSFrechetDistance (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double* dist);

/**
* Calculate the
* [Frechet distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance)
* between two geometries,
* a similarity measure for linear features. For more precision, first
* densify the inputs.
* \param[in] g1 Input geometry
* \param[in] g2 Input geometry
* \param[in] densifyFrac The largest % of the overall line length that
*            any given two-point segment should be
* \param[out] dist Pointer to be filled in with distance result
* \return 1 on success, 0 on exception.
* \see geos::algorithm::distance::DiscreteFrechetDistance
*/
int GEOSFrechetDistanceDensify (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2,
    double densifyFrac,
    double* dist);

/**
* Calculate the length of a LineString.
* Only works for LineString inputs, returns exception otherwise.
*
* \param[in] g Input geometry
* \param[out] length Pointer to be filled in with length result
* \return 1 on success, 0 on exception.
*/
int GEOSGeomGetLength (const(GEOSGeometry)* g, double* length);

/**
* The closest points of the two geometries.
* The first point comes from g1 geometry and the second point comes from g2.
*
* \param[in] g1 Input geometry
* \param[in] g2 Input geometry
* \return A coordinate sequence with the two points, or NULL on exception.
* Caller must free with GEOSCoordSeq_destroy().
*/
GEOSCoordSequence* GEOSNearestPoints (
    const(GEOSGeometry)* g1,
    const(GEOSGeometry)* g2);

/* ========== Algorithms ========== */

/**
* For the points formed by the six input ordinates,
* walking from A to B and then to P.
* \param Ax X coordinate of A
* \param Ay Y coordinate of A
* \param Bx X coordinate of B
* \param By Y coordinate of B
* \param Px X coordinate of P
* \param Py Y coordinate of P
* \return  -1 if reaching P takes a counter-clockwise (left) turn,
*           1 if reaching P takes a clockwise (right) turn,
*           0 if P is collinear with A-B
*/
int GEOSOrientationIndex (
    double Ax,
    double Ay,
    double Bx,
    double By,
    double Px,
    double Py);

/* ========= Reader and Writer APIs ========= */

/* ========= WKT Reader ========= */

/**
* Allocate a new \ref GEOSWKTReader.
* \returns a new reader. Caller must free with GEOSWKTReader_destroy()
*/
GEOSWKTReader* GEOSWKTReader_create ();

/**
* Free the memory associated with a \ref GEOSWKTReader.
* \param reader The reader to destroy.
*/
void GEOSWKTReader_destroy (GEOSWKTReader* reader);

/**
* Use a reader to parse the well-known text representation of
* a geometry, and return an allocated geometry.
* \param reader A WKT reader object, caller retains ownership
* \param wkt The WKT string to parse, caller retains ownership
* \return A \ref GEOSGeometry, caller to free with GEOSGeom_destroy())
*/
GEOSGeometry* GEOSWKTReader_read (GEOSWKTReader* reader, const(char)* wkt);

/* ========= WKT Writer ========= */

/**
* Allocate a new \ref GEOSWKTWriter.
* \returns a new writer. Caller must free with GEOSWKTWriter_destroy()
*/
GEOSWKTWriter* GEOSWKTWriter_create ();

/**
* Free the memory associated with a \ref GEOSWKTWriter.
* \param writer The writer to destroy.
*/
void GEOSWKTWriter_destroy (GEOSWKTWriter* writer);

/**
* Writes out the well-known text representation of a geometry,
* using the trim, rounding and dimension settings of the writer.
* \param writer A \ref GEOSWKTWriter.
* \param g Input geometry
* \return A newly allocated string containing the WKT output or NULL on exception.
* Caller must free with GEOSFree()
*/
char* GEOSWKTWriter_write (GEOSWKTWriter* writer, const(GEOSGeometry)* g);

/**
* Sets the number trimming option on a \ref GEOSWKTWriter.
* With trim set to 1, the writer will strip trailing 0's from
* the output coordinates. With 0, all coordinates will be
* padded with 0's out to the rounding precision.
* \param writer A \ref GEOSWKTWriter.
* \param trim The trimming behaviour to set, 1 for 'on', 0 for 'off', default 'off'
*/
void GEOSWKTWriter_setTrim (GEOSWKTWriter* writer, char trim);

/**
* Sets the number places after the decimal to output in
* WKT.
* \param writer A \ref GEOSWKTWriter.
* \param precision The desired precision, default 16.
*/
void GEOSWKTWriter_setRoundingPrecision (GEOSWKTWriter* writer, int precision);

/**
* Sets whether or not to write out XY or XYZ coordinates.
* Legal values are 2 or 3.
* \param writer A \ref GEOSWKTWriter.
* \param dim The desired dimension, default 2.
*/
void GEOSWKTWriter_setOutputDimension (GEOSWKTWriter* writer, int dim);

/**
* Reads the current output dimension from a \ref GEOSWKTWriter.
* \param writer A \ref GEOSWKTWriter.
* \return The current dimension.
*/
int GEOSWKTWriter_getOutputDimension (GEOSWKTWriter* writer);

/**
* Sets the format for 3D outputs. The "old 3D" format does not
* include a dimensionality tag, eg. "POINT(1 2 3)" while the new (ISO)
* format does includes a tag, eg "POINT Z (1 2 3)".
* \param writer A \ref GEOSWKTWriter.
* \param useOld3D True to use the old format, false is the default.
*/
void GEOSWKTWriter_setOld3D (GEOSWKTWriter* writer, int useOld3D);

/* ========== WKB Reader ========== */

/**
* Allocate a new \ref GEOSWKBReader.
* \returns a new reader. Caller must free with GEOSWKBReader_destroy()
*/
GEOSWKBReader* GEOSWKBReader_create ();

/**
* Free the memory associated with a \ref GEOSWKBReader.
* \param reader The reader to destroy.
*/
void GEOSWKBReader_destroy (GEOSWKBReader* reader);

/**
* Read a geometry from a well-known binary buffer.
* \param reader A \ref GEOSWKBReader
* \param wkb A pointer to the buffer to read from
* \param size The number of bytes of data in the buffer
* \return A \ref GEOSGeometry built from the WKB, or NULL on exception.
*/
GEOSGeometry* GEOSWKBReader_read (
    GEOSWKBReader* reader,
    const(ubyte)* wkb,
    size_t size);

/**
* Read a geometry from a **hex encoded** well-known binary buffer.
* \param reader A \ref GEOSWKBReader
* \param hex A pointer to the buffer to read from
* \param size The number of bytes of data in the buffer
* \return A \ref GEOSGeometry built from the HEX WKB, or NULL on exception.
*/
GEOSGeometry* GEOSWKBReader_readHEX (
    GEOSWKBReader* reader,
    const(ubyte)* hex,
    size_t size);

/* ========== WKB Writer ========== */

/**
* Allocate a new \ref GEOSWKBWriter.
* \returns a new writer. Caller must free with GEOSWKBWriter_destroy()
*/
GEOSWKBWriter* GEOSWKBWriter_create ();

/**
* Free the memory associated with a \ref GEOSWKBWriter.
* \param writer The writer to destroy.
*/
void GEOSWKBWriter_destroy (GEOSWKBWriter* writer);

/**
* Write out the WKB representation of a geometry.
* \param writer The \ref GEOSWKBWriter controlling the
* writing.
* \param g Geometry to convert to WKB
* \param size Pointer to write the size of the final output WKB to
* \return The WKB representation. Caller must free with GEOSFree()
*/
ubyte* GEOSWKBWriter_write (
    GEOSWKBWriter* writer,
    const(GEOSGeometry)* g,
    size_t* size);

/**
* Write out the **hex** WKB representation of a geometry.
* \param writer The \ref GEOSWKBWriter controlling the
* writing.
* \param g Geometry to convert to WKB
* \param size Pointer to write the size of the final output WKB to
* \return The HEX WKB representation. Caller must free with GEOSFree()
*/
ubyte* GEOSWKBWriter_writeHEX (
    GEOSWKBWriter* writer,
    const(GEOSGeometry)* g,
    size_t* size);

/**
* Read the current output dimension of the writer.
* Either 2 or 3 dimensions.
* Return current number of dimensions.
* \param writer The writer to read from.
* \return Number of dimensions (2 or 3)
*/
int GEOSWKBWriter_getOutputDimension (const(GEOSWKBWriter)* writer);

/**
* Set the output dimensionality of the writer. Either
* 2 or 3 dimensions.
* \param writer The writer to read from.
* \param newDimension The dimensionality desired
*/
void GEOSWKBWriter_setOutputDimension (GEOSWKBWriter* writer, int newDimension);

/**
* Find whether the writer will use WKB
* [byte order](https://en.wikipedia.org/wiki/Endianness)
* that is big or little endian.
* The return value is a member of \ref GEOSWKBByteOrders.
* \param writer The writer to read byte order from
* \return The current byte order
*/
int GEOSWKBWriter_getByteOrder (const(GEOSWKBWriter)* writer);

/**
* Set the output byte order of the writer, using
* a value from \ref GEOSWKBByteOrders enum.
* \param writer The writer to set byte order on
* \param byteOrder Desired byte order
*/
void GEOSWKBWriter_setByteOrder (GEOSWKBWriter* writer, int byteOrder);

/**
* Find whether the writer will use
* [WKB](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry#Well-known_binary)
* that is ISO flavor or "extended" flavor. The flavor
* determines how extra dimensionality is encoded with the
* type number, and whether SRID can be included in the WKB.
* ISO flavor does not support SRID embedding. ISO flavor
* is "more standard" for 3D output. GEOS can read both flavors.
* The return value is a member of \ref GEOSWKBFlavors.
* \param writer The writer to read flavor from
* \return The current flavor
*/
int GEOSWKBWriter_getFlavor (const(GEOSWKBWriter)* writer);

/**
* Set the output flavor of the writer, using
* a value from \ref GEOSWKBFlavors enum.
* \param writer The writer to set flavor on
* \param flavor Desired flavor
*/
void GEOSWKBWriter_setFlavor (GEOSWKBWriter* writer, int flavor);

/**
* Read the current SRID embedding value from the writer.
* \param writer The writer to check SRID value on
*/
char GEOSWKBWriter_getIncludeSRID (const(GEOSWKBWriter)* writer);

/**
* Specify whether SRID values should be output in WKB.
* Many WKB readers do not support SRID values, use with caution.
* \param writer The writer to set SRID output on
* \param writeSRID Set to 1 to include SRID, 0 otherwise
*/
void GEOSWKBWriter_setIncludeSRID (GEOSWKBWriter* writer, const char writeSRID);

/**
* Free strings and byte buffers returned by functions such
* as GEOSWKBWriter_write(),
* GEOSWKBWriter_writeHEX() and GEOSWKTWriter_write(), etc.
* \param buffer The memory to free
*/
void GEOSFree (void* buffer);

/* ========= GeoJSON Reader ========= */

/**
* Allocate a new \ref GEOSGeoJSONReader.
* \returns a new reader. Caller must free with GEOSGeoJSONReader_destroy()
*/
GEOSGeoJSONReader* GEOSGeoJSONReader_create ();

/**
* Free the memory associated with a \ref GEOSGeoJSONReader.
* \param reader The reader to destroy.
*/
void GEOSGeoJSONReader_destroy (GEOSGeoJSONReader* reader);

/**
* Use a reader to parse a GeoJSON. A single geometry or feature is
* converted into a geometry. A featurecollection is converted into a
* geometrycollection. Feature properties are not read.
* \param reader A GeoJSON reader object, caller retains ownership
* \param geojson The json string to parse, caller retains ownership
* \return A \ref GEOSGeometry, caller to free with GEOSGeom_destroy())
*/
GEOSGeometry* GEOSGeoJSONReader_readGeometry (
    GEOSGeoJSONReader* reader,
    const(char)* geojson);

/* ========= GeoJSON Writer ========= */

/**
* Allocate a new \ref GEOSGeoJSONWriter.
* \returns a new writer. Caller must free with GEOSGeoJSONWriter_destroy()
*/
GEOSGeoJSONWriter* GEOSGeoJSONWriter_create ();

/**
* Free the memory associated with a \ref GEOSGeoJSONWriter.
* \param writer The writer to destroy.
*/
void GEOSGeoJSONWriter_destroy (GEOSGeoJSONWriter* writer);

/**
* Write out the GeoJSON representation of a geometry. Note that writing a GeoJSON
* Feature or FeatureCollection is unsupported through the GEOS C API.
* \param writer A GeoJSON reader object, caller retains ownership.
* \param g The geometry to convert, caller retains ownership.
* \param indent The indentation used. Use -1 for no formatting.
* \return A char pointer, caller to free with GEOSFree())
*/
char* GEOSGeoJSONWriter_writeGeometry (
    GEOSGeoJSONWriter* writer,
    const(GEOSGeometry)* g,
    int indent);

/* #ifndef GEOS_USE_ONLY_R_API */

/* ====================================================================== */
/* DEPRECATIONS */
/* ====================================================================== */

/**
* \deprecated in 3.3.0, use GEOSOffsetCurve() instead
*/
GEOSGeometry* GEOSSingleSidedBuffer (
    const(GEOSGeometry)* g,
    double width,
    int quadsegs,
    int joinStyle,
    double mitreLimit,
    int leftSide);

/**
* \deprecated in 3.3.0, use GEOSOffsetCurve() instead
*/
GEOSGeometry* GEOSSingleSidedBuffer_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    double width,
    int quadsegs,
    int joinStyle,
    double mitreLimit,
    int leftSide);

/**
* \deprecated in 3.5.0. Use GEOS_init_r() and set the message handlers using
* GEOSContext_setNoticeHandler_r() and/or GEOSContext_setErrorHandler_r()
*/
GEOSContextHandle_t initGEOS_r (
    GEOSMessageHandler notice_function,
    GEOSMessageHandler error_function);

/**
* \deprecated in 3.5.0, replaced by GEOS_finish_r()
*/
void finishGEOS_r (GEOSContextHandle_t handle);

/**
* \deprecated use \ref GEOSWKTReader and GEOSWKTReader_read_r()
*/
GEOSGeometry* GEOSGeomFromWKT_r (GEOSContextHandle_t handle, const(char)* wkt);

/**
* \deprecated use \ref GEOSWKTWriter and GEOSWKTWriter_write_r()
*/
char* GEOSGeomToWKT_r (GEOSContextHandle_t handle, const(GEOSGeometry)* g);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_getOutputDimension_r()
*/
int GEOS_getWKBOutputDims_r (GEOSContextHandle_t handle);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_setOutputDimension_r()
*/
int GEOS_setWKBOutputDims_r (GEOSContextHandle_t handle, int newDims);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_getByteOrder_r()
*/
int GEOS_getWKBByteOrder_r (GEOSContextHandle_t handle);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_setByteOrder_r()
*/
int GEOS_setWKBByteOrder_r (GEOSContextHandle_t handle, int byteOrder);

/**
* \deprecated use \ref GEOSWKBReader and GEOSWKBReader_read_r()
*/
GEOSGeometry* GEOSGeomFromWKB_buf_r (
    GEOSContextHandle_t handle,
    const(ubyte)* wkb,
    size_t size);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_write_r()
*/
ubyte* GEOSGeomToWKB_buf_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    size_t* size);

/**
* \deprecated use \ref GEOSWKBReader and GEOSWKBReader_readHEX_r()
*/
GEOSGeometry* GEOSGeomFromHEX_buf_r (
    GEOSContextHandle_t handle,
    const(ubyte)* hex,
    size_t size);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_writeHEX_r()
*/
ubyte* GEOSGeomToHEX_buf_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g,
    size_t* size);

/**
* \deprecated use \ref GEOSWKTReader and GEOSWKTReader_read_r()
*/
GEOSGeometry* GEOSGeomFromWKT (const(char)* wkt);

/**
* \deprecated use \ref GEOSWKTWriter and GEOSWKTWriter_write()
*/
char* GEOSGeomToWKT (const(GEOSGeometry)* g);

/**
* \deprecated use \ref GEOSWKBWriter and GEOS_getWKBOutputDims()
*/
int GEOS_getWKBOutputDims ();

/**
* \deprecated use \ref GEOSWKBWriter and GEOS_setWKBOutputDims()
*/
int GEOS_setWKBOutputDims (int newDims);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_getByteOrder()
*/
int GEOS_getWKBByteOrder ();

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_setByteOrder()
*/
int GEOS_setWKBByteOrder (int byteOrder);

/**
* \deprecated use \ref GEOSWKBReader and GEOSWKBWriter_read()
*/
GEOSGeometry* GEOSGeomFromWKB_buf (const(ubyte)* wkb, size_t size);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_write()
*/
ubyte* GEOSGeomToWKB_buf (const(GEOSGeometry)* g, size_t* size);

/**
* \deprecated use \ref GEOSWKBReader and GEOSWKBWriter_readHEX()
*/
GEOSGeometry* GEOSGeomFromHEX_buf (const(ubyte)* hex, size_t size);

/**
* \deprecated use \ref GEOSWKBWriter and GEOSWKBWriter_writeHEX()
*/
ubyte* GEOSGeomToHEX_buf (const(GEOSGeometry)* g, size_t* size);

/**
* \deprecated in 3.3.0: use GEOSUnaryUnion() instead
*/
GEOSGeometry* GEOSUnionCascaded (const(GEOSGeometry)* g);

/**
* \deprecated in 3.3.0: use GEOSUnaryUnion_r() instead
*/
GEOSGeometry* GEOSUnionCascaded_r (
    GEOSContextHandle_t handle,
    const(GEOSGeometry)* g);

/* ====================================================================== */
/* END DEPRECATIONS */
/* ====================================================================== */

// extern "C"

/* #ifndef GEOS_C_H_INCLUDED */
