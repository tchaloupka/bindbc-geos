/**
 * More D idiomatic API wrapper around GEOS C API
 */
module bindbc.geos.wrapper;

import core.stdc.stdarg;
import core.stdc.stdio;
import core.stdc.stdlib;
import core.stdc.string : strlen;
import bindbc.geos.libgeos;
import std.algorithm : among;
import std.math : isNaN;

enum MessageLevel { notice, error }

alias MessageHandler = void function(MessageLevel, const(char)[]) nothrow;

/**
 * Initializes GEOS for the current thread.
 * It must be called for each tred prior to the usage.
 * If no message handler is provided, default'll be used that just outputs messages to the STDOUT/STDERR.
 *
 * `finishGEOS` must be called to free resources on thread exit
 */
void initGEOS(MessageHandler handler = null) @trusted nothrow @nogc
{
    assert(ctx is null, "GEOS context already initialized for this thread");
    ctx = GEOS_init_r();
    if (ctx is null) assert(0, "Failed to initialize GEOS");
    if (handler !is null) msgHandler = handler;
    GEOSContext_setNoticeHandler_r(ctx, &msgHandlerImpl!(MessageLevel.notice));
    GEOSContext_setErrorHandler_r(ctx, &msgHandlerImpl!(MessageLevel.error));
    wktWriter = GEOSWKTWriter_create_r(ctx);
    if (wktWriter is null) assert(0, "Failed to create WKT writer");
    GEOSWKTWriter_setTrim_r(ctx, wktWriter, 1); // enable trailing 0 trimming
}

/// Cleans up the current GEOS context
void finishGEOS() @trusted nothrow @nogc
{
    assert(ctx !is null, "GEOS context not initialized");
    assert(wktWriter !is null, "WKT writer not initialized");
    GEOSWKTWriter_destroy_r(ctx, wktWriter);
    GEOS_finish_r(ctx);
}

/// Returns the current GEOS version string. eg: "3.10.0"
string geosVersion() @trusted nothrow @nogc
{
    import std.exception : assumeUnique;
    auto pc = GEOSversion();
    return pc[0..strlen(pc)].assumeUnique;
}

version (unittest)
{
    static this() @trusted nothrow @nogc
    {
        initGEOS();
        printf("Initialized GEOS Version: %s\n", geosVersion().ptr); // points to a zero terminated string literal
    }
    static ~this() @safe nothrow @nogc { finishGEOS(); }
}

/// Simple type used to work with point coordinates
union PointXYZ(int dims) if (dims.among(2, 3))
{
    /// Point constructor
    this(double[dims] coords)
    {
        this.coords = coords;
    }

    static if (dims == 3)
    {
        /// Ditto
        this(double x, double y, double z = double.init)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }
    else
    {
        /// Ditto
        this(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
    }

    double[dims] coords;
    struct {
        double x, y;
        static if (dims == 3) double z;
    }
}

alias PointZ = PointXYZ!3;
alias Point = PointXYZ!2;

/// Geometry extent coordinates (rectangle)
struct Extent
{
    double xmin;
    double ymin;
    double xmax;
    double ymax;
}

/// Simple wrapper for GEOS allocated strings.
struct GeosString
{
    @disable this(this); /// disable copying

    ~this() nothrow @nogc @trusted
    {
        if (str) GEOSFree_r(ctx, c);
    }

    const(char)[] str() const nothrow @nogc @trusted
    {
        if (c is null) return null;
        return c[0..len];
    }

    private:
    char* c;
    size_t len;

    this(char* str) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (str !is null, "Null pointer provided")
    {
        this.c = str;
        this.len = strlen(str);
    }
}

/**
 * Wrapper around GEOS library `GEOSCoordSequence` type.
 * Underlying instance is automatically freed using RAII.
 */
struct CoordSequence
{
// int 	GEOSCoordSeq_copyToBuffer (const GEOSCoordSequence *s, double *buf, int hasZ, int hasM)
// int 	GEOSCoordSeq_copyToArrays (const GEOSCoordSequence *s, double *x, double *y, double *z, double *m)

    /**
     * Constructs `CoordSequence` with known size but with uninitialized coordinates.
     *
     * Note: Coordinate values of this sequence'll be left undefined and must be set separately.
     * It'll contain random garbage values otherwise. Thats why it's marked as `@system`.
     */
    this(uint size, uint dims) @system nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (size > 0, "Can't create empty sequence")
    in (dims == 2 || dims == 3, "Dimensionality of the coordinates must be in 2 or 3")
    {
        seq = GEOSCoordSeq_create_r(ctx, size, dims);
        if (seq is null) assert(0, "Failed to initialize CoordSequence");
    }

    /// Constructs `CoordSequence` from individual array for each dimension.
    this(const(double)[] x, const(double)[] y, const(double)[] z = null) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (x.length && y.length, "x and y arrays must have elements")
    in (x.length == y.length && (!z.length || x.length == z.length), "Coordinate arrays must be of equal length")
    in (x.length <= uint.max, "Sequence too large")
    {
        static if (geosSupport >= GEOSSupport.geos_3_10)
        {
            seq = GEOSCoordSeq_copyFromArrays_r(ctx, x.ptr, y.ptr, z.ptr, null, cast(uint)x.length);
            if (seq is null) assert(0, "Failed to initialize CoordSequence");
        }
        else
        {
            if (!z.length)
            {
                this(cast(uint)x.length, 2);
                for (uint i=0; i<x.length; ++i) {
                    this[i, 0] = x[i];
                    this[i, 1] = y[i];
                }
            } else  {
                this(cast(uint)x.length, 3);
                for (uint i=0; i<x.length; ++i) {
                    this[i, 0] = x[i];
                    this[i, 1] = y[i];
                    this[i, 2] = z[i];
                }
            }
        }
    }

    /// Constructs `CoordSequence` from array of supported point types
    this(P)(const(P)[] points) @trusted nothrow @nogc
        if (isPoint!P)
    in (ctx !is null, "GEOS thread context not initialized")
    in (points.length <= uint.max, "Too many points")
    in (points.length > 0, "Can't create empty sequence")
    {
        static if (isStaticArrayPoint!P) {
            seq = GEOSCoordSeq_copyFromBuffer_r(ctx, points.ptr, points.length, points[0].length-2, 0);
            if (seq is null) assert(0, "Failed to initialize CoordSequence");
        }
        else static if (isDynamicArrayPoint!P)
        {
            // get required dimensions
            auto dims = points[0].length;
            assert(dims.among(2, 3), "Invalid coordinate dimensions");

            // prepare space for coordinates
            seq = GEOSCoordSeq_create_r(ctx, cast(uint)points.length, cast(uint)dims);
            if (seq is null) assert(0, "Failed to initialize CoordSequence");

            // assign individual points
            foreach (i, ref p; points) this[i] = p;
        }
        else
        {
            // prepare space for coordinates
            static if (is(typeof(points[0].z) : double)) enum dims = 3;
            else enum dims = 2;
            seq = GEOSCoordSeq_create_r(ctx, cast(uint)points.length, dims);
            if (seq is null) assert(0, "Failed to initialize CoordSequence");

            // assign individual points
            foreach (i, ref p; points) this[i] = p;
        }
    }

    /// ditto - just variadic variant for coordinate structs as a source
    this(P)(const(P)[] points...) @trusted nothrow @nogc
        if (isPointStruct!P)
    {
        this(points);
    }

    @disable this(this); /// Copy is not allowed - use move() or clone()

    ~this() @trusted nothrow @nogc
    {
        if (seq !is null && own)
            GEOSCoordSeq_destroy_r(ctx, seq);
    }

    /// Create a new copy of the geometry.
    CoordSequence clone() const @trusted nothrow @nogc
    {
        CoordSequence cseq;
        if (this.seq) {
            cseq.seq = GEOSCoordSeq_clone_r(ctx, this.seq);
            if (cseq.seq is null) assert(0, "Failed to clone sequence");
        }
        return cseq;
    }

    /// Copy the contents of a coordinate sequence to individual dimension buffers
    void toArrays(double[] x, double[] y, double[] z = null) const @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    in (x.length == y.length && (z.length == 0 || z.length == x.length) && x.length == length(), "Lengths don't match")
    {
        if (z.length)
        {
            assert(dimensions() == 3, "Incompatible dimensions");
            static if (geosSupport >= GEOSSupport.geos_3_10) {
                auto res = GEOSCoordSeq_copyToArrays_r(ctx, seq, x.ptr, y.ptr, z.ptr, null);
                assert(res == 1, "Failed to copy coordinates");
            }
            else
            {
                for (uint i=0; i<length(); ++i) {
                    auto res = GEOSCoordSeq_getXYZ_r(ctx, seq, i, &x[i], &y[i], &z[i]);
                    assert(res == 1, "Failed to copy coordinates");
                }
            }
            return;
        }

        static if (geosSupport >= GEOSSupport.geos_3_10) {
            auto res = GEOSCoordSeq_copyToArrays_r(ctx, seq, x.ptr, y.ptr, null, null);
            assert(res == 1, "Failed to copy coordinates");
        }
        else
        {
            for (uint i=0; i<length(); ++i) {
                auto res = GEOSCoordSeq_getXY_r(ctx, seq, i, &x[i], &y[i]);
                assert(res == 1, "Failed to copy coordinates");
            }
        }
    }

    /// Copy the contents of a coordinate sequence to a buffer of doubles (XYXY or XYZXYZ)
    void toBuffer(int N)(double[N][] buf) const @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    in (buf.length == length(), "Lengths don't match")
    {
        static if (N == 2)
        {
            static if (geosSupport >= GEOSSupport.geos_3_10) {
                auto res = GEOSCoordSeq_copyToBuffer_r(ctx, seq, cast(double*)buf.ptr, 0, 0);
                assert(res == 1, "Failed to copy coordinates");
            } else {
                for (uint i=0; i<length(); ++i) {
                    auto res = GEOSCoordSeq_getXY_r(ctx, seq, i, &buf[i][0], &buf[i][1]);
                    assert(res == 1, "Failed to copy coordinates");
                }
            }
        }
        else
        {
            assert(dimensions() == 3, "Incompatible dimensions");
            static if (geosSupport >= GEOSSupport.geos_3_10) {
                auto res = GEOSCoordSeq_copyToBuffer_r(ctx, seq, cast(double*)buf.ptr, 1, 0);
                assert(res == 1, "Failed to copy coordinates");
            } else {
                for (uint i=0; i<length(); ++i) {
                    auto res = GEOSCoordSeq_getXYZ_r(ctx, seq, i, &buf[i][0], &buf[i][1], &buf[i][2]);
                    assert(res == 1, "Failed to copy coordinates");
                }
            }
        }
    }

    /// Get size info from a coordinate sequence.
    uint length() const @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    {
        uint sz = void;
        auto r = GEOSCoordSeq_getSize_r(ctx, seq, &sz);
        assert(r == 1, "Failed to get sequence size");
        return sz;
    }

    /// Get dimension info from a coordinate sequence.
    uint dimensions() const @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    {
        uint dims = void;
        auto r = GEOSCoordSeq_getDimensions_r(ctx, seq, &dims);
        assert(r == 1, "Failed to get sequence dimensions");
        return dims;
    }

    /**
     * Check orientation of a coordinate sequence. Closure of the sequence is assumed. Invalid (collapsed) sequences will return false.
     * Sequence must have at least 4 points.
     */
    bool isCCW() const @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    in (length >= 4, "not enough points to check CCW")
    {
        char isccw = void;
        auto r = GEOSCoordSeq_isCCW_r(ctx, seq, &isccw);
        assert(r == 1, "Failed to check if sequence is counter clockwise");
        return isccw == 1;
    }

    /// Gets value of requested coordinate in sequence
    double opIndex(size_t idx, size_t dim) const @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    in (idx < length, "invalid index")
    in (dim < dimensions, "invalid dimension")
    {
        double v = void;
        auto r = GEOSCoordSeq_getOrdinate_r(ctx, seq, cast(uint)idx, cast(uint)dim, &v);
        assert (r == 1, "Failed to get coordinate value");
        return v;
    }

    /// ditto
    const(PointZ) opIndex(size_t idx) const @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    in (idx < length, "invalid index")
    {
        PointZ p = void;
        int r;
        if (dimensions == 2) {
            r = GEOSCoordSeq_getXY_r(ctx, seq, cast(uint)idx, &p.x, &p.y);
            p.z = double.init;
        }
        else r = GEOSCoordSeq_getXYZ_r(ctx, seq, cast(uint)idx, &p.x, &p.y, &p.z);
        assert (r == 1, "Failed to get point coordinates");
        return p;
    }

    /// Support for dollar indexing
    int opDollar(size_t pos)() const @trusted nothrow @nogc
    {
        static if (pos == 0) return length();
        else static if (pos == 1) return dimensions();
        else static assert(0, "Invalid index dimension");
    }

    /// Set Nth ordinate value in a coordinate sequence.
    double opIndexAssign(double value, size_t idx, size_t dim) @trusted nothrow @nogc
    in (seq !is null, "unitialized sequence")
    in (idx < length, "invalid index")
    in (dim < dimensions, "invalid dimension")
    {
        auto r = GEOSCoordSeq_setOrdinate_r(ctx, seq, cast(uint)idx, cast(uint)dim, value);
        assert (r == 1, "Failed to set coordinate value");
        return value;
    }

    /// ditto for generic acceptable point types
    auto opIndexAssign(P)(auto ref P point, size_t idx) @trusted nothrow @nogc
        if (isPoint!P)
    in (seq !is null, "unitialized sequence")
    in (idx < length, "invalid index")
    {
        int res;
        static if (isStaticArrayPoint!P)
        {
            assert(dimensions() == point.length, "Incompatible dimensions");
            static if (N == 2)
                res = GEOSCoordSeq_setXY_r(ctx, seq, cast(uint)idx, cast(double)point[0], cast(double)point[1]);
            else
                res = GEOSCoordSeq_setXYZ_r(ctx, seq, cast(uint)idx, cast(double)point[0], cast(double)point[1], cast(double)point[2]);
        }
        else static if (isDynamicArrayPoint!P)
        {
            assert(point.length.among(2, 3), "Invalid coordinate length");
            assert(dimensions() == point.length, "Invalid coordinate length");
            if (dimensions() == 2)
                res = GEOSCoordSeq_setXY_r(ctx, seq, cast(uint)idx, cast(double)point[0], cast(double)point[1]);
            else
                res = GEOSCoordSeq_setXYZ_r(ctx, seq, cast(uint)idx, cast(double)point[0], cast(double)point[1], cast(double)point[2]);
        }
        else
        {
            static if (is(typeof(point.z) : double)) {
                assert(dimensions() == 3, "Incompatible dimensions");
                res = GEOSCoordSeq_setXYZ_r(ctx, seq, cast(uint)idx, cast(double)point.x, cast(double)point.y, cast(double)point.z);
            } else {
                assert(dimensions() == 2, "Incompatible dimensions");
                res = GEOSCoordSeq_setXY_r(ctx, seq, cast(uint)idx, cast(double)point.x, cast(double)point.y);
            }
        }
        assert (res == 1, "Failed to set coordinate value");
        return point;
    }

    /// Returns RAW GEOS handle to it's CoordSequence type
    GEOSCoordSequence* handle() @safe nothrow @nogc
    {
        return seq;
    }

    private:
    GEOSCoordSequence* seq;
    bool own = true;
}

///
@("CoordSequence")
@safe unittest
{
    // construct with separate dimension arrays
    auto coords = CoordSequence([-1, 1, 1, -1, -1], [-2, -2, 2, 2, -2]);
    assert(coords.dimensions == 2);
    assert(coords.length == 5);
    assert(coords.clone().length == 5);
    assert(coords.isCCW);
    assert(coords[0, 0] == -1);
    assert(coords[0, 1] == -2);
    assert(coords[$-1, 0] == -1);
    assert(coords[$-1, $-1] == -2);
    assert(coords[1].x == 1);
    assert(coords[1].z.isNaN);
    coords[0, 0] = 42;
    assert(coords[0, 0] == 42);

    // construct with uninitialized coordinates
    () @trusted
    {
        coords = CoordSequence(1, 2);
        coords[0] = [42, 666];
        assert(coords[0].x == 42);
        assert(coords[0].y == 666);
    }();

    // construct with point array
    coords = CoordSequence([[1,2], [3,4]]);
    assert(coords[0].x == 1);
    assert(coords[0].y == 2);
    assert(coords[1].x == 3);
    assert(coords[1].y == 4);

    coords = CoordSequence([Point(1,2), Point(3,4)]);
    assert(coords[0].x == 1);
    assert(coords[0].y == 2);
    assert(coords[1].x == 3);
    assert(coords[1].y == 4);

    // copy to buffers
    double[] x = new double[2];
    double[] y = new double[2];
    double[2][] p = new double[2][2];
    coords.toArrays(x, y);
    assert(x == [1, 3]);
    assert(y == [2, 4]);
    coords.toBuffer(p);
    assert(p == [[1, 2], [3, 4]]);
}

/**
 * Wrapper around GEOS library `GEOSGeometry*` type.
 * Underlying instance is automatically freed using RAII.
 */
struct Geometry
{
    /**
     * Gets type of the geometry from `GEOSGeomTypes` enum.
     * Geometry must be already initialized.
     * Uses: `GEOSGeomTypeId` function
     */
    GEOSGeomTypes typeId() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto tid = GEOSGeomTypeId_r(ctx, g);
        assert(tid >= 0, "Error getting geometry typeId");
        return cast(GEOSGeomTypes)tid;
    }

    /// Create a new copy of the geometry.
    Geometry clone() const @trusted nothrow @nogc
    {
        Geometry c;
        if (this.g) {
            c.g = GEOSGeom_clone_r(ctx, this.g);
            if (c.g is null) assert(0, "Failed to clone geometry");
        }
        return c;
    }

    alias x = getPointCoordinate!GEOSGeomGetX_r; /// Returns the X coordinate of a Point geometry
    alias y = getPointCoordinate!GEOSGeomGetY_r; /// Returns the Y coordinate of a Point geometry
    alias z = getPointCoordinate!GEOSGeomGetZ_r; /// Returns the Z coordinate of a Point geometry

    alias xmin = getPointCoordinate!GEOSGeom_getXMin_r; /// Finds the minimum X value in the geometry.
    alias xmax = getPointCoordinate!GEOSGeom_getXMax_r; /// Finds the maximum X value in the geometry.
    alias ymin = getPointCoordinate!GEOSGeom_getYMin_r; /// Finds the minimum Y value in the geometry.
    alias ymax = getPointCoordinate!GEOSGeom_getYMax_r; /// Finds the maximum Y value in the geometry.

    /**
     * Return the planar dimensionality of the geometry.
     *   0 for point, multipoint
     *   1 for linestring, multilinestring
     *   2 for polygon, multipolygon
     */
    int dimensions() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        return GEOSGeom_getDimensions_r(ctx, g);
    }

    /// Finds the extent (minimum and maximum X and Y value) of the geometry.
    void extent(out double xmin, out double ymin, out double xmax, out double ymax) const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        static if (geosSupport >= GEOSSupport.geos_3_11)
        {
            auto r = GEOSGeom_getExtent_r(ctx, g, &xmin, &ymin, &xmax, &ymax);
            if (r != 1) assert(0, "Error getting extent for geometry");
        }
        else
        {
            xmin = this.xmin();
            ymin = this.ymin();
            xmax = this.xmax();
            ymax = this.ymax();
        }
    }

    /// ditto
    Extent extent() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        Extent ext;
        extent(ext.xmin, ext.ymin, ext.xmax, ext.ymax);
        return ext;
    }

    /// Get the total number of points in a geometry, of any type.
    uint coordNum() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto r = GEOSGetNumCoordinates_r(ctx, g);
        assert(r != -1, "Failed to get number of coordinates");
        return cast(uint)r;
    }

    /// Return the coordinate sequence underlying the given geometry (Must be a LineString,
    /// LinearRing or Point).
    CoordSequence coordSeq() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        CoordSequence cseq;
        cseq.seq = cast(GEOSCoordSequence*)GEOSGeom_getCoordSeq_r(ctx, g);
        assert(cseq.seq !is null, "Failed to get coordinate sequence for geometry");
        cseq.own = false;
        return cseq;
    }

    /**
     * Return the cartesian dimension of the geometry.
     *
     *   * 2 for XY data
     *   * 3 for XYZ data
     */
    uint coordDimensions() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto ret = GEOSGeom_getCoordinateDimension_r(ctx, g);
        assert(ret.among(2,3), "Failed to determine geometry dimensions");
        return cast(uint)ret;
    }

    @disable this(this); /// Copy is not allowed - use move() or clone()

    ~this() @trusted nothrow @nogc
    {
        if (g && own) GEOSGeom_destroy_r(ctx, g);
    }

    /// Creates an empty Point geometry.
    static Geometry createEmptyPoint() @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    {
        Geometry res;
        res.g = GEOSGeom_createEmptyPoint_r(ctx);
        if (res.g is null) assert(0, "Failed to create empty Point");
        return res;
    }

    /// Creates a Point geometry from a pair of coordinates.
    static Geometry createPoint(double x, double y) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    {
        Geometry res;
        res.g = GEOSGeom_createPointFromXY_r(ctx, x, y);
        if (res.g is null) assert(0, "Failed to create Point");
        return res;
    }

    /// Creates a Point geometry from `CoordSequence`
    static Geometry createPoint(CoordSequence coord) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (coord.seq, "Uninitialized coord sequence")
    {
        Geometry res;
        res.g = GEOSGeom_createPoint_r(ctx, coord.seq);
        if (res.g is null) assert(0, "Failed to create Point");
        coord.seq = null; // geometry takes ownership of coord sequence
        return res;
    }

    /// Creates LinearRing geometry from `CoordSequence`
    static Geometry createLinearRing(CoordSequence coords) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (coords.seq, "Uninitialized coord sequence")
    {
        Geometry res;
        res.g = GEOSGeom_createLinearRing_r(ctx, coords.seq);
        if (res.g is null) assert(0, "Failed to create LinearRing");
        coords.seq = null; // geometry takes ownership of coords sequence
        return res;
    }

    /// Creates LineString geometry from `CoordSequence`
    static Geometry createLineString(CoordSequence coords) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (coords.seq, "Uninitialized coord sequence")
    {
        Geometry res;
        res.g = GEOSGeom_createLineString_r(ctx, coords.seq);
        if (res.g is null) assert(0, "Failed to create LineString");
        coords.seq = null; // geometry takes ownership of coords sequence
        return res;
    }

    /// Creates an empty LineString.
    static Geometry createEmptyLineString() @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    {
        Geometry res;
        res.g = GEOSGeom_createEmptyLineString_r(ctx);
        if (res.g is null) assert(0, "Failed to create empty LineString");
        return res;
    }

    /// Creates an empty Polygon.
    static Geometry createEmptyPolygon() @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    {
        Geometry res;
        res.g = GEOSGeom_createEmptyPolygon_r(ctx);
        if (res.g is null) assert(0, "Failed to create empty Polygon");
        return res;
    }

    /// Creates a Polygon geometry from line ring geometries.
    /// Created Geometry takes ownership of the supplied ones.
    static Geometry createPolygon(Geometry shell, Geometry[] holes...) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (shell.g, "Uninitialized shell geometry")
    {
        import core.exception : onOutOfMemoryError;
        Geometry res;
        if (holes.length)
        {
            // we need to make an array of Geometry pointers first
            if (holes.length >= 64)
            {
                GEOSGeometry*[64] buf;
                foreach (i, ref h; holes) {
                    assert(h.g, "Uninitialized hole geometry");
                    buf[i] = h.g;
                }
                res.g = GEOSGeom_createPolygon_r(ctx, shell.g, buf.ptr, cast(uint)holes.length);
            }
            else
            {
                auto pholes = cast(GEOSGeometry**)malloc((GEOSGeometry*).sizeof * holes.length);
                if (!pholes) onOutOfMemoryError();
                foreach (i, ref h; holes) {
                    assert(h.g, "Uninitialized hole geometry");
                    pholes[i] = h.g;
                }
                res.g = GEOSGeom_createPolygon_r(ctx, shell.g, pholes, cast(uint)holes.length);
                free(pholes);
            }
        }
        else res.g = GEOSGeom_createPolygon_r(ctx, shell.g, null, 0);
        if (res.g is null) assert(0, "Failed to create Polygon");
        shell.g = null; // geometry takes ownership
        foreach (ref h; holes) h.g = null; // for holes too
        return res;
    }

    /**
     * Create a geometry collection.
     * Created Geometry takes ownership of the supplied ones.
     */
    static Geometry createCollection(GEOSGeomTypes type, Geometry[] geoms...) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (geoms.length, "No source geometry provided")
    {
        import core.exception : onOutOfMemoryError;
        Geometry res;

        // we need to make an array of Geometry pointers first
        if (geoms.length >= 64)
        {
            GEOSGeometry*[64] buf;
            foreach (i, ref h; geoms) {
                assert(h.g, "Uninitialized geometry");
                buf[i] = h.g;
            }
            res.g = GEOSGeom_createCollection_r(ctx, type, buf.ptr, cast(uint)geoms.length);
        }
        else
        {
            auto pgeoms = cast(GEOSGeometry**)malloc((GEOSGeometry*).sizeof * geoms.length);
            if (!pgeoms) onOutOfMemoryError();
            foreach (i, ref h; geoms) {
                assert(h.g, "Uninitialized geometry");
                pgeoms[i] = h.g;
            }
            res.g = GEOSGeom_createCollection_r(ctx, type, pgeoms, cast(uint)geoms.length);
            free(pgeoms);
        }
        if (res.g is null) assert(0, "Failed to create geometry collection");
        foreach (ref h; geoms) h.g = null; // resulting geometry collection gets ownership of source ones
        return res;
    }

    /// Creates an empty geometry collection.
    static Geometry createEmptyCollection(GEOSGeomTypes type) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    {
        Geometry res;
        res.g = GEOSGeom_createEmptyCollection_r(ctx, type);
        if (res.g is null) assert(0, "Failed to create empty geometry collection");
        return res;
    }

    /// Create a rectangular polygon from bounding coordinates. Will return a point geometry if width and height are 0.
    static Geometry createRectangle(double xmin, double ymin, double xmax, double ymax) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (!xmin.isNaN && !ymin.isNaN && !xmax.isNaN && !ymax.isNaN, "Undefined coordinate provided")
    {
        static if (geosSupport >= GEOSSupport.geos_3_11) {
            Geometry res;
            res.g = GEOSGeom_createRectangle_r(ctx, type);
            if (res.g is null) assert(0, "Failed to create empty geometry collection");
            return res;
        }
        else {
            if (xmin == xmax && ymin == ymax)
                return createPoint(xmin, ymin);
            else
                return createPolygon(createLinearRing(CoordSequence(
                    Point(xmin, ymin),
                    Point(xmax, ymin),
                    Point(xmax, ymax),
                    Point(xmin, ymax),
                    Point(xmin, ymin)
                )));
        }
    }

    /// Return the anonymous "user data" for this geometry. User data must be managed by the caller,
    /// and freed before the geometry is freed.
    void* userData() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        return GEOSGeom_getUserData_r(ctx, g);
    }

    /// Set the anonymous "user data" for this geometry. Don't forget to free the user data before
    /// freeing the geometry.
    void userData(void* data) @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        GEOSGeom_setUserData_r(ctx, g, data);
    }

    /// Returns the number of sub-geometries immediately under a multi-geometry or collection or 1
    /// for a simple geometry. For nested collections, remember to check if returned sub-geometries
    /// are themselves also collections.
    int numGeometries() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        return GEOSGetNumGeometries_r(ctx, g);
    }

    /// Returns the specified sub-geometry of a collection. For a simple geometry, returns itself.
    /// Warning: returned geometry must not be used after original geometry is discarded (as internal pointer would be invalid)
    Geometry geometry(int idx) const @trusted nothrow @nogc
    in (idx < numGeometries(), "Invalid index")
    {
        Geometry ret;
        ret.g = cast(GEOSGeom_t*)GEOSGetGeometryN_r(ctx, g, idx);
        assert(ret.g !is null, "Failed to get indexed subgeometry");
        ret.own = false;
        return ret;
    }

    /**
     * Organize the elements, rings, and coordinate order of geometries in a consistent way, so that
     * geometries that represent the same object can be easily compared. Modifies the geometry
     * in-place.
     *
     * Normalization ensures the following:
     *
     *   * Lines are oriented to have smallest coordinate first (apart from duplicate endpoints)
     *   * Rings start with their smallest coordinate (using XY ordering)
     *   * Polygon shell rings are oriented CW, and holes CCW
     *   * Collection elements are sorted by their first coordinate
     *
     * Use before calling `equalsExact` to avoid false "not equal" results.
     */
    void normalize() @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto r = GEOSNormalize_r(ctx, g);
        assert(r == 0, "Failed to normalize geometry");
    }

    static if (geosSupport >= GEOSSupport.geos_3_11)
    {
        /**
        * Works from start of each coordinate sequence in the geometry, retaining points that are
        * further away from the previous retained point than the tolerance value.
        *
        * Removing repeated points with a non-zero tolerance may result in an invalid geometry being
        * returned. Be sure to check and repair validity.
        *
        * Returns:
        *     A geometry with all points within the tolerance of each other removed.
        */
        Geometry removeRepeatedpoints(double tolerance = 0.0) const @trusted nothrow @nogc
        in (g !is null, "unitialized geometry")
        {
            Geometry ret;
            ret.g = GEOSRemoveRepeatedPoints_r(ctx, g, tolerance);
            assert(ret.g !is null, "Failed to remove repeated points");
            return ret;
        }
    }

    /// Tests whether the input geometry is empty. If the geometry or any component is non-empty,
    /// the geometry is non-empty. An empty geometry has no boundary or interior.
    bool isEmpty() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto ret = GEOSisEmpty_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is empty");
        return ret == 1;
    }

    /// Tests whether the input geometry is a ring. Rings are linestrings, without
    /// self-intersections, with start and end point being identical.
    bool isRing() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto ret = GEOSisRing_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is ring");
        return ret == 1;
    }

    /// Tests whether the input geometry has z coordinates.
    bool hasZ() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto ret = GEOSHasZ_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry has Z");
        return ret == 1;
    }

    /// Tests whether the input geometry is closed. A closed geometry is a linestring or
    /// multilinestring with the start and end points being the same.
    bool isClosed() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto ret = GEOSisClosed_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is closed");
        return ret == 1;
    }

    /// Tests whether the input geometry is "simple". Mostly relevant for linestrings. A "simple"
    /// linestring has no self-intersections.
    bool isSimple() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto ret = GEOSisSimple_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is simple");
        return ret == 1;
    }

    /**
     * Check the validity of the provided geometry.
     *
     *   * All points are valid.
     *   * All non-zero-length linestrings are valid.
     *   * Polygon rings must be non-self-intersecting, and interior rings contained within exterior rings.
     *   * Multi-polygon components may not touch or overlap.
     */
    bool isValid() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto ret = GEOSisValid_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is valid");
        return ret == 1;
    }

    /// Return the human readable reason a geometry is invalid, "Valid Geometry" string otherwise.
    GeosString isValidReason() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        auto c = GEOSisValidReason_r(ctx, g);
        assert(c !is null, "Failed to check if geometry is valid");
        return GeosString(c);
    }

    /// In one step, calculate and return the validity, the human readable validity reason and a
    /// point at which validity rules are broken.
    bool isValidDetail(out GeosString reason, out Geometry location) const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        char* c;
        GEOSGeometry* g;
        auto ret = GEOSisValidDetail_r(ctx, g, 0, &c, &g);
        assert(ret.among(0, 1), "Failed to check if geometry is valid");
        if (c) reason = GeosString(c);
        if (g) location.g = g;
        return ret == 1;
    }

    /// True if geometries cover the same space on the place.
    bool equals(const Geometry other) const @trusted nothrow @nogc
    in (g !is null && other.g !is null, "unitialized geometry")
    {
        auto r = GEOSEquals_r(ctx, this.g, other.g);
        assert(r.among(0, 1), "Failed to determine if geometries are equal");
        return r == 1;
    }

    /// Determine pointwise equivalence of two geometries, by checking that they have identical
    /// structure and that each vertex of g2 is within the distance tolerance of the corresponding
    /// vertex in g1. Unlike GEOSEquals(), geometries that are topologically equivalent but have
    /// different representations (e.g., LINESTRING (0 0, 1 1) and MULTILINESTRING ((0 0, 1 1)) )
    /// are not considered equal by GEOSEqualsExact().
    bool equalsExact(const Geometry other, double tolerance = 0.0) const @trusted nothrow @nogc
    in (g !is null && other.g !is null, "unitialized geometry")
    {
        auto r = GEOSEqualsExact_r(ctx, this.g, other.g, tolerance);
        assert(r.among(0, 1), "Failed to determine if geometries are equal");
        return r == 1;
    }

    /// Generates geometry as [WKT](https://libgeos.org/specifications/wkt/) string
    void toString(S)(auto ref S sink) const @trusted
    in (g !is null, "unitialized geometry")
    {
        char* s = GEOSWKTWriter_write_r(ctx, wktWriter, g);
        if (s is null) assert(0, "Failed to generate WKT from geometry");
        sink(s[0..strlen(s)]);
        GEOSFree_r(ctx, s);
    }

    /// ditto
    string toString() const @trusted
    in (g !is null, "unitialized geometry")
    {
        char* s = GEOSWKTWriter_write_r(ctx, wktWriter, g);
        if (s is null) assert(0, "Failed to generate WKT from geometry");
        string str = s[0..strlen(s)].idup;
        GEOSFree_r(ctx, s);
        return str;
    }

    /// Returns RAW GEOS handle to it's Geometry type
    GEOSGeometry* handle() @safe nothrow @nogc
    {
        return g;
    }

    /// Calculate the area of a geometry.
    double area() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        double area;
        auto ret = GEOSArea_r(ctx, g, &area);
        assert(ret == 1, "Error calculating geometry area");
        return area;
    }

    /// Calculate the length of a geometry.
    double length() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        double len;
        auto ret = GEOSLength_r(ctx, g, &len);
        assert(ret == 1, "Error calculating geometry length");
        return len;
    }

    private:
    GEOSGeometry* g;
    bool own = true;

    double getPointCoordinate(alias fn)() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    in (typeId() == GEOSGeomTypes.GEOS_POINT, "Operation is valid only on Point geometry")
    {
        double coord = void;
        auto r = fn(ctx, g, &coord);
        assert(r == 1, "Error getting coordinate");
        return coord;
    }
}

///
@("Point geometry")
@safe
unittest
{
    import core.lifetime : move;
    import std.array : Appender;

    // Point from x, y
    auto point = Geometry.createPoint(1, 2);
    assert(point.g !is null);
    assert(point.typeId == GEOSGeomTypes.GEOS_POINT);
    auto s = point.toString();
    assert(s == "POINT (1 2)", s);

    assert(point.x == 1);
    assert(point.y == 2);
    assert(point.z.isNaN);
    assert(point.coordNum == 1);
    assert(point.dimensions == 0);
    assert(point.xmin == 1);
    assert(point.xmax == 1);
    assert(point.ymin == 2);
    assert(point.ymax == 2);

    auto ext = point.extent();
    assert(ext.xmin == 1);
    assert(ext.xmax == 1);
    assert(ext.ymin == 2);
    assert(ext.ymax == 2);

    () @trusted
    {
        int data = 42;
        point.userData = &data;
        assert (point.userData == &data);
    }();

    // WKT string generation
    Appender!string buf;
    point.clone().toString((const(char)[] s) {buf ~= s; });
    assert(buf.data == "POINT (1 2)", buf.data);

    // Point from coord sequence
    auto cseq = CoordSequence(Point(1,2));
    point = Geometry.createPoint(cseq.move); // move is needed as CoordSequence is not copyable and Geometry takes ownership of the sequence
    assert(cseq.seq is null);
    assert(point.x == 1);
    assert(point.y == 2);
    assert(!point.isEmpty);

    // Empty point
    point = Geometry.createEmptyPoint();
    assert(point.isEmpty);
}

///
@("Linestring geometry")
@safe
unittest
{
    auto g = Geometry.createLineString(CoordSequence(Point(1,2), Point(3,4)));
    assert(g.typeId == GEOSGeomTypes.GEOS_LINESTRING);
    // assert(g[0].x == 1);
    // assert(g[0].y == 2);
    // assert(g[1].x == 3);
    // assert(g[1].y == 4);
}

/// Corresponds with PostGIS output numbers
enum LineCrossingDirection
{
    noCross = 0,
    left = -1,
    right = 1,
    multicrossLeft = -2,
    multicrossRight = 2,
    multicrossFirstLeft = -3,
    multicrossFirstRight = 3
}

/**
 * Determines direction in which the line A is crossed by line B.
 * See: [ST_LineCrossingDirection](https://postgis.net/docs/ST_LineCrossingDirection.html)
 */
LineCrossingDirection lineCrossingDirection(ref const Geometry lineA, ref const Geometry lineB)
{
    //TODO: https://github.com/postgis/postgis/blob/f6def67654c25d812446239036cee44812613748/liblwgeom/lwalgorithm.c
    // https://github.com/postgis/postgis/blob/45e491242d9d63beb7eac0961b16439d00c93bd4/liblwgeom/cunit/cu_algorithm.c
    return LineCrossingDirection.noCross;
}

private:
GEOSContextHandle_t ctx; // thread GEOS context
GEOSWKTWriter* wktWriter;
MessageHandler msgHandler = &defaultMessageHandler;

void defaultMessageHandler(MessageLevel lvl, const(char)[] msg) @trusted nothrow @nogc
{
    final switch (lvl)
    {
        case MessageLevel.notice:
            fprintf(stdout, "%s\n", msg.ptr);
            break;
        case MessageLevel.error:
            fprintf(stderr, "%s\n", msg.ptr);
            break;
    }
}

extern (C) nothrow @trusted
void msgHandlerImpl(MessageLevel lvl)(const(char)* fmt, ...)
{
    va_list ap;

    // /* Determine required size. */
    // va_start(ap, fmt);
    // auto n = vsnprintf(null, 0, fmt, ap);
    // va_end(ap);

    // if (n < 0) return;

    // size_t size = n + 1;      /* One extra byte for '\0' */
    // auto p = cast(char*)malloc(size);
    // if (p is null) return;

    // va_start(ap, fmt);
    // n = vsnprintf(p, size, fmt, ap);
    // va_end(ap);

    // if (n < 0) {
    //     free(p);
    //     return;
    // }

    // Workaround: abowe doesn't work due to the https://issues.dlang.org/show_bug.cgi?id=21425
    char[512] buf;
    va_start(ap, fmt);
    auto n = vsnprintf(buf.ptr, buf.length, fmt, ap);
    va_end(ap);

    msgHandler(lvl, buf[0..n]);
    // free(p);
}

template isStaticArrayPoint(P)
{
    enum isStaticArrayPoint = is(P : C[N],C,N) && is(C : double) (N == 2 || N == 3);
}

template isDynamicArrayPoint(P)
{
    enum isDynamicArrayPoint = !isPointStruct!P && is(P : C[], C) && is(C : double);
}

template isPointStruct(P)
{
    enum isPointStruct = is(typeof(P.x) : double) && is(typeof(P.y) : double);
}

template isPoint(P)
{
    enum isPoint = isStaticArrayPoint!P || isDynamicArrayPoint!P || isPointStruct!P;
}

static assert(isPoint!Point);
static assert(!isStaticArrayPoint!Point);
static assert(!isDynamicArrayPoint!Point);
