/**
 * More D idiomatic API wrapper around GEOS C API
 */
module bindbc.geos.wrapper;

import core.lifetime : move;
import core.stdc.stdarg;
import core.stdc.stdio;
import core.stdc.stdlib;
import core.stdc.string : strlen;
import bindbc.geos.libgeos;
import std.algorithm : among, max, min;
import core.exception : onOutOfMemoryError;
import std.math : isNaN;
import std.meta : allSatisfy;
import std.traits : allSameType, isNumeric;
public import bindbc.geos.libgeos : GEOSGeomTypes;

enum MessageLevel { notice, error }

alias MessageHandler = void function(MessageLevel, const(char)[]) nothrow;

/**
 * Initializes GEOS for the current thread.
 * It must be called for each thread prior to the usage.
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
    wktReader = GEOSWKTReader_create_r(ctx);
    if (wktReader is null) assert(0, "Failed to create WKT reader");
    wktWriter = GEOSWKTWriter_create_r(ctx);
    if (wktWriter is null) assert(0, "Failed to create WKT writer");
    GEOSWKTWriter_setTrim_r(ctx, wktWriter, 1); // enable trailing 0 trimming
}

/// Cleans up the current GEOS context
void finishGEOS() @trusted nothrow @nogc
{
    assert(ctx !is null, "GEOS context not initialized");
    assert(wktReader !is null, "WKT reader not initialized");
    assert(wktWriter !is null, "WKT writer not initialized");
    GEOSWKTReader_destroy_r(ctx, wktReader);
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

        /// Converts to 2D point
        Point point2d() @safe pure nothrow @nogc
        {
            return Point(x, y);
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

    auto opAssign(P)(auto ref P p) if (is(P : PointZ) || is(P : Point))
    {
        this.x = p.x;
        this.y = p.y;
        static if (dims == 3) {
            static if (is(P : PointZ))
                this.z = p.z;
            else
                this.z = double.init;
        }
        return this;
    }

    double[dims] coords;
    struct {
        double x, y;
        static if (dims == 3) double z;
    }

    T opCast(T)() const if (is(T == Geometry))
    {
        return Geometry.createPoint(this);
    }
}

alias PointZ = PointXYZ!3;
alias Point = PointXYZ!2;

@("Point")
@safe unittest
{
    auto pt = Point(0.1, 0.2);
    auto geom = cast(Geometry)pt;
    assert(geom.coordDimensions == 2);
    assert(geom.x == pt.x);
    assert(geom.y == pt.y);
}

/// Geometry extent coordinates (rectangle)
struct Extent
{
    double xmin;
    double ymin;
    double xmax;
    double ymax;

    T opCast(T)() const if (is(T == Geometry))
    {
        return Geometry.createRectangle(xmin, ymin, xmax, ymax);
    }
}

@("Extent")
@safe unittest
{
    auto ext = Extent(0, 0, 1, 1);
    assert((cast(Geometry)ext).extent == ext);
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

    /// ditto - just variadic variant for coordinate structures as a source
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
    in (seq !is null, "uninitialized sequence")
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
    in (seq !is null, "uninitialized sequence")
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
    in (seq !is null, "uninitialized sequence")
    {
        uint sz = void;
        auto r = GEOSCoordSeq_getSize_r(ctx, seq, &sz);
        assert(r == 1, "Failed to get sequence size");
        return sz;
    }

    /// Get dimension info from a coordinate sequence.
    uint dimensions() const @trusted nothrow @nogc
    in (seq !is null, "uninitialized sequence")
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
    in (seq !is null, "uninitialized sequence")
    in (length >= 4, "not enough points to check CCW")
    {
        char isccw = void;
        auto r = GEOSCoordSeq_isCCW_r(ctx, seq, &isccw);
        assert(r == 1, "Failed to check if sequence is counter clockwise");
        return isccw == 1;
    }

    /// Gets value of requested coordinate in sequence
    double opIndex(size_t idx, size_t dim) const @trusted nothrow @nogc
    in (seq !is null, "uninitialized sequence")
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
    in (seq !is null, "uninitialized sequence")
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
    in (seq !is null, "uninitialized sequence")
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
    in (seq !is null, "uninitialized sequence")
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
                if (dimensions() == 2)
                    res = GEOSCoordSeq_setXY_r(ctx, seq, cast(uint)idx, cast(double)point.x, cast(double)point.y);
                else
                    res = GEOSCoordSeq_setXYZ_r(ctx, seq, cast(uint)idx, cast(double)point.x, cast(double)point.y, cast(double)point.z);
            } else {
                assert(dimensions() == 2, "Incompatible dimensions");
                res = GEOSCoordSeq_setXY_r(ctx, seq, cast(uint)idx, cast(double)point.x, cast(double)point.y);
            }
        }
        assert (res == 1, "Failed to set coordinate value");
        return point;
    }

    /// Returns new CoordSequence instance from sliced coordinates
    CoordSequence opSlice(uint start, uint end) @trusted nothrow @nogc
    in (seq !is null, "uninitialized sequence")
    in (start >= 0 && end <= length, "invalid index")
    {
        CoordSequence res = CoordSequence(end - start, this.dimensions);
        for (size_t i = start; i<end; i++)
            res[i-start] = this[i];
        return res;
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

    // slice
    auto sliced = coords[1..3];
    assert(sliced.length == 2);
    assert(sliced.dimensions == 2);
    assert(sliced[0, 0] == 1);
    assert(sliced[0, 1] == -2);
    assert(sliced[1, 0] == 1);
    assert(sliced[1, 1] == 2);

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
    in (g !is null, "uninitialized geometry")
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

    alias x = getCoordinate!GEOSGeomGetX_r; /// Returns the X coordinate of a Point geometry
    alias y = getCoordinate!GEOSGeomGetY_r; /// Returns the Y coordinate of a Point geometry
    alias z = getCoordinate!GEOSGeomGetZ_r; /// Returns the Z coordinate of a Point geometry

    alias xmin = getCoordinate!GEOSGeom_getXMin_r; /// Finds the minimum X value in the geometry.
    alias xmax = getCoordinate!GEOSGeom_getXMax_r; /// Finds the maximum X value in the geometry.
    alias ymin = getCoordinate!GEOSGeom_getYMin_r; /// Finds the minimum Y value in the geometry.
    alias ymax = getCoordinate!GEOSGeom_getYMax_r; /// Finds the maximum Y value in the geometry.

    /**
     * Return the planar dimensionality of the geometry.
     *   0 for point, multipoint
     *   1 for linestring, multilinestring
     *   2 for polygon, multipolygon
     */
    int dimensions() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        return GEOSGeom_getDimensions_r(ctx, g);
    }

    /// Finds the extent (minimum and maximum X and Y value) of the geometry.
    void extent(out double xmin, out double ymin, out double xmax, out double ymax) const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
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
    in (g !is null, "uninitialized geometry")
    {
        Extent ext;
        extent(ext.xmin, ext.ymin, ext.xmax, ext.ymax);
        return ext;
    }

    /// Get the total number of points in a geometry, of any type.
    uint coordNum() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        auto r = GEOSGetNumCoordinates_r(ctx, g);
        assert(r != -1, "Failed to get number of coordinates");
        return cast(uint)r;
    }

    /// Return the coordinate sequence underlying the given geometry (Must be a LineString,
    /// LinearRing or Point).
    CoordSequence coordSeq() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
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
    in (g !is null, "uninitialized geometry")
    {
        auto ret = GEOSGeom_getCoordinateDimension_r(ctx, g);
        assert(ret.among(2,3), "Failed to determine geometry dimensions");
        return cast(uint)ret;
    }

    /**
     * Initializes geometry from Well Known Text (WKT) format.
     * If creation fails, resulting Geometry would be uninitialized.
     */
    this(const(char)[] wkt) nothrow @nogc @trusted
    in (ctx !is null, "GEOS thread context not initialized")
    in (wkt.length, "Empty WKT")
    {
        import std.internal.cstring : tempCString;
        this.g = GEOSWKTReader_read_r(ctx, wktReader, wkt.tempCString);
    }

    @disable this(this); /// Copy is not allowed - use move() or clone()

    ~this() @trusted nothrow @nogc
    {
        if (g && own) GEOSGeom_destroy_r(ctx, g);
    }

    /// Checks if the Geometry is initialized
    bool opCast(T)() const if (is(T == bool))
    {
        return g !is null;
    }

    // Sets point coordinates
    ref Geometry opAssign(P)(auto ref P pt) @trusted
        if (is(P : PointZ) || is(P : Point))
    {
        if (g is null) this = Geometry.createPoint(pt);
        else
        {
            assert(this.typeId == GEOSGeomTypes.GEOS_POINT, "No point geometry");
            auto cs = GEOSGeom_getCoordSeq_r(ctx, g);
            int r;
            static if (is(P == Point))
                r = GEOSCoordSeq_setXY_r(ctx, cast(GEOSCoordSequence*)cs, 0, pt.x, pt.y);
            else
            {
                if (dimensions == 2)
                    r = GEOSCoordSeq_setXY_r(ctx, cast(GEOSCoordSequence*)cs, 0, pt.x, pt.y);
                else
                    r = GEOSCoordSeq_setXYZ_r(ctx, cast(GEOSCoordSequence*)cs, 0, pt.x, pt.y, pt.z);
            }
            assert(r != 0, "GEOSCoordSeq_setXY()");
        }
        return this;
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
    static Geometry createPoint(double x, double y, double z = double.init) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    in (!x.isNaN && !y.isNaN, "Invalid coordinates")
    {
        Geometry res;
        if (z.isNaN)
        {
            res.g = GEOSGeom_createPointFromXY_r(ctx, x, y);
            if (res.g is null) assert(0, "Failed to create Point");
            return res.move();
        }
        else
        {
            auto coords = CoordSequence(1, 3);
            coords[0] = PointZ(x, y, z);
            return createPoint(coords.move);
        }
    }

    /// Creates a Point geometry from a pair of coordinates.
    static Geometry createPoint(P)(auto ref P pt) @trusted nothrow @nogc
        if (isPoint!P)
    {
        static if (isPointStruct!P)
        {
            static if (is(typeof(pt.z) : double)) return createPoint(pt.x, pt.y, pt.z);
            else return createPoint(pt.x, pt.y);
        }
        else
        {
            if (pt.length == 3) return createPoint(pt[0], pt[1], pt[2]);
            else return createPoint(pt[0], pt[1]);
        }
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

    /// Create LinearRing geometry from points
    static Geometry createLinearRing(P)(const(P)[] points) @trusted if (isPoint!P)
    {
        return Geometry.createLinearRing(CoordSequence(points));
    }

    /// ditto
    static Geometry createLinearRing(P)(const(P)[] points...) if (isPointStruct!P)
    {
        return Geometry.createLinearRing(points);
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

    /// Creates LineString geometry from points
    static Geometry createLineString(P)(const(P)[] points) @trusted if (isPoint!P)
    {
        return Geometry.createLineString(CoordSequence(points));
    }

    /// ditto
    static Geometry createLineString(P)(const(P)[] points...) if (isPointStruct!P)
    {
        return Geometry.createLineString(points);
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
    static Geometry createPolygon(ARGS...)(Geometry shell, auto ref ARGS holes) @trusted nothrow @nogc
        if (ARGS.length == 0 || (allSameType!ARGS && is(ARGS[0] == Geometry)))
    in (ctx !is null, "GEOS thread context not initialized")
    in (shell.g, "Uninitialized shell geometry")
    {
        Geometry res;
        static if (ARGS.length)
        {
            // we need to make an array of Geometry pointers first
            static if (holes.length >= 64)
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

    /// Creates polygon from points
    static Geometry createPolygon(P)(const(P)[] points) @trusted if (isPoint!P)
    {
        auto g = Geometry.createPolygon(Geometry.createLinearRing(CoordSequence(points)));
        if (!g.isValid) {
            auto vg = GEOSMakeValid_r(ctx, g.g);
            GEOSGeom_destroy_r(ctx, g.g);
            g.g = vg;
        }
        return g;
    }

    /// ditto
    static Geometry createPolygon(P)(const(P)[] points...) if (isPointStruct!P)
    {
        return Geometry.createPolygon(points);
    }

    /**
     * Create a geometry collection.
     * Created Geometry takes ownership of the supplied ones.
     */
    static Geometry createCollection(ARGS...)(GEOSGeomTypes type, auto ref ARGS geoms) @trusted nothrow @nogc
        if (ARGS.length >= 1 && allSameType!ARGS && is(ARGS[0] == Geometry))
    in (ctx !is null, "GEOS thread context not initialized")
    in (geoms.length, "No source geometry provided")
    {
        Geometry res;

        // we need to make an array of Geometry pointers first
        static if (ARGS.length <= 64)
        {
            GEOSGeometry*[ARGS.length] buf;
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
    in (g !is null, "uninitialized geometry")
    {
        return GEOSGeom_getUserData_r(ctx, g);
    }

    /// Set the anonymous "user data" for this geometry. Don't forget to free the user data before
    /// freeing the geometry.
    void userData(void* data) @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        GEOSGeom_setUserData_r(ctx, g, data);
    }

    /// Returns the number of sub-geometries immediately under a multi-geometry or collection or 1
    /// for a simple geometry. For nested collections, remember to check if returned sub-geometries
    /// are themselves also collections.
    int numGeometries() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
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
    in (g !is null, "uninitialized geometry")
    {
        auto r = GEOSNormalize_r(ctx, g);
        assert(r == 0, "Failed to normalize geometry");
    }

    /// Returns the union of all components of a single geometry. Usually used to convert a
    /// collection into the smallest set of polygons that cover the same area.
    Geometry union_() @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        Geometry ret;
        ret.g = GEOSUnaryUnion_r(ctx, g);
        assert(ret.g, "GEOSUnaryUnion()");
        return ret;
    }

    /// Returns newly allocated geometry of the difference
    Geometry difference()(auto ref const(Geometry) geom) @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    in (geom.g !is null, "uninitialized geometry")
    {
        Geometry ret;
        ret.g = GEOSDifference_r(ctx, g, geom.g);
        assert(ret.g, "GEOSDifference()");
        return ret;
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
        Geometry removeRepeatedPoints(double tolerance = 0.0) const @trusted nothrow @nogc
        in (g !is null, "uninitialized geometry")
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
    in (g !is null, "uninitialized geometry")
    {
        auto ret = GEOSisEmpty_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is empty");
        return ret == 1;
    }

    /// Tests whether the input geometry is a ring. Rings are linestrings, without
    /// self-intersections, with start and end point being identical.
    bool isRing() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        auto ret = GEOSisRing_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is ring");
        return ret == 1;
    }

    /// Tests whether the input geometry has z coordinates.
    bool hasZ() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        auto ret = GEOSHasZ_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry has Z");
        return ret == 1;
    }

    /// Tests whether the input geometry is closed. A closed geometry is a linestring or
    /// multilinestring with the start and end points being the same.
    bool isClosed() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        auto ret = GEOSisClosed_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is closed");
        return ret == 1;
    }

    /// Tests whether the input geometry is "simple". Mostly relevant for linestrings. A "simple"
    /// linestring has no self-intersections.
    bool isSimple() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
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
    in (g !is null, "uninitialized geometry")
    {
        auto ret = GEOSisValid_r(ctx, g);
        assert(ret.among(0, 1), "Failed to check if geometry is valid");
        return ret == 1;
    }

    /// Return the human readable reason a geometry is invalid, "Valid Geometry" string otherwise.
    GeosString isValidReason() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        auto c = GEOSisValidReason_r(ctx, g);
        assert(c !is null, "Failed to check if geometry is valid");
        return GeosString(c);
    }

    /// In one step, calculate and return the validity, the human readable validity reason and a
    /// point at which validity rules are broken.
    bool isValidDetail(out GeosString reason, out Geometry location) const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
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
    in (g !is null && other.g !is null, "uninitialized geometry")
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
    in (g !is null && other.g !is null, "uninitialized geometry")
    {
        auto r = GEOSEqualsExact_r(ctx, this.g, other.g, tolerance);
        assert(r.among(0, 1), "Failed to determine if geometries are equal");
        return r == 1;
    }

    /// Generates geometry as [WKT](https://libgeos.org/specifications/wkt/) string
    void toString(S)(auto ref S sink) const @trusted
    in (g !is null, "uninitialized geometry")
    {
        import std.range : put;
        char* s = GEOSWKTWriter_write_r(ctx, wktWriter, g);
        if (s is null) assert(0, "Failed to generate WKT from geometry");
        put(sink, s[0..strlen(s)]);
        GEOSFree_r(ctx, s);
    }

    /// ditto
    string toString() const @trusted
    in (g !is null, "uninitialized geometry")
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
    in (g !is null, "uninitialized geometry")
    {
        double area;
        auto ret = GEOSArea_r(ctx, g, &area);
        assert(ret == 1, "Error calculating geometry area");
        return area;
    }

    /// Calculate the length of a geometry.
    double length() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
    {
        double len;
        auto ret = GEOSLength_r(ctx, g, &len);
        assert(ret == 1, "Error calculating geometry length");
        return len;
    }

    /**
     * Appends point or linestring to the Geometry.
     * If geometry is uninitialized, it would be initialized with the provided geometry.
     * Otherwise provided point(s) would be appended to the end of the current ones.
     *
     * Only valid for Point or Linestring geometries.
     */
    void opOpAssign(string op, T)(auto ref T value) @trusted if (op=="~" && (isPoint!T || is(T == Geometry)))
    in (!this || typeId.among(GEOSGeomTypes.GEOS_POINT, GEOSGeomTypes.GEOS_LINESTRING), "Unsupported operation")
    {
        static if (is(T == Geometry)) {
            assert(
                value && value.typeId.among(GEOSGeomTypes.GEOS_POINT, GEOSGeomTypes.GEOS_LINESTRING),
                "Invalid input geometry");
        }

        if (!this)
        {
            static if (is(T == Geometry)) this = value.clone();
            else this = createPoint(value);
        }
        else
        {
            static if (isPoint!T)
            {
                auto newCoords = CoordSequence(this.coordNum+1, this.coordDimensions);
                auto coords = this.coordSeq;
                for (uint i=0; i<coords.length; i++)
                    newCoords[i] = coords[i];
                newCoords[coords.length] = value;
                this = Geometry.createLineString(newCoords.move());
            }
            else
            {
                auto newCoords = CoordSequence(this.coordNum+value.coordNum, this.coordDimensions);
                auto coords = this.coordSeq;
                auto gcoords = value.coordSeq;
                for (uint i; i<coords.length; i++)
                    newCoords[i] = coords[i];
                for (uint i; i<gcoords.length; i++)
                    newCoords[coords.length + i] = gcoords[i];
                this = Geometry.createLineString(newCoords.move);
            }
        }
    }

    /// True if geometry is completely within geom, and not touching the boundary of geom
    bool within(ref Geometry geom) const @trusted nothrow @nogc
    {
        auto r = GEOSWithin_r(ctx, this.g, geom.g);
        assert(r != 2, "GEOSWithin()");
        return r == 1;
    }

    /// True if geom is completely within geometry.
    bool contains(ref Geometry geom) const @trusted nothrow @nogc
    {
        auto r = GEOSContains_r(ctx, this.g, geom.g);
        assert(r != 2, "GEOSContains()");
        return r == 1;
    }

    /// Calculates distance between two geometries
    double distance(ref Geometry geom) const @trusted nothrow @nogc
    {
        double dist;
        auto r = GEOSDistance_r(ctx, this.g, geom.g, &dist);
        assert(r == 1, "GEOSDistance()");
        return dist;
    }

    private:
    GEOSGeometry* g;
    bool own = true;

    double getCoordinate(alias fn)() const @trusted nothrow @nogc
    in (g !is null, "uninitialized geometry")
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
    assert(g.xmin == 1);
    assert(g.ymin == 2);
    assert(g.xmax == 3);
    assert(g.ymax == 4);
}

@("Linestring append")
@safe
unittest
{
    Geometry a;
    a ~= Point(1,2);
    auto s = a.toString();
    assert(s == "POINT (1 2)", s);
    a ~= Point(3, 4);
    s = a.toString();
    assert(s == "LINESTRING (1 2, 3 4)", s);

    Geometry b = Geometry("LINESTRING(11 22, 33 44)");
    a ~= b;
    s = a.toString();
    assert(s == "LINESTRING (1 2, 3 4, 11 22, 33 44)", s);
    s = b.toString();
    assert(s == "LINESTRING (11 22, 33 44)", s);
    a = Geometry.createPoint(1, 2);
    a ~= b;
    s = a.toString();
    assert(s == "LINESTRING (1 2, 11 22, 33 44)", s);

    a = Geometry.init;
    a ~= b;
    s = a.toString();
    assert(s == "LINESTRING (11 22, 33 44)", s);
}

@("Geometry - within/contains")
@safe unittest
{
    auto pt = Geometry.createPoint(0.5, 0.5);
    auto poly = Geometry.createPolygon(Geometry.createLinearRing(CoordSequence(
        Point(0, 0),
        Point(1, 0),
        Point(1, 1),
        Point(0, 1),
        Point(0, 0)
    )));

    assert(poly.contains(pt));
    assert(!pt.contains(poly));
    assert(pt.within(poly));
    assert(!poly.within(pt));
}

@("Geometry - polygon")
@safe unittest
{
    auto shell = Geometry.createLinearRing(
        Point(0, 0),
        Point(1, 0),
        Point(1, 1),
        Point(0, 1),
        Point(0, 0)
    );

    auto hole = Geometry.createLinearRing(
        Point(0.4, 0.4),
        Point(0.6, 0.4),
        Point(0.6, 0.6),
        Point(0.4, 0.6),
        Point(0.4, 0.4)
    );

    auto poly = Geometry.createPolygon(shell.clone(), hole.clone());
    auto s = poly.toString();
    assert(s == "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0), (0.4 0.4, 0.6 0.4, 0.6 0.6, 0.4 0.6, 0.4 0.4))", s);

    poly = Geometry.createPolygon(shell.clone());
    auto phole = Geometry.createPolygon(hole.clone());
    poly = poly.difference(phole); // another way to make a hole
    assert(phole); // not consumed
    s = poly.toString();
    assert(s == "POLYGON ((0 1, 1 1, 1 0, 0 0, 0 1), (0.6 0.4, 0.6 0.6, 0.4 0.6, 0.4 0.4, 0.6 0.4))", s);

    poly = Geometry.createPolygon(shell.clone());
    poly = poly.difference(Geometry.createPolygon(
        Point(2, 0),
        Point(3, 0),
        Point(3, 1),
        Point(2, 1),
        Point(2, 0)
    ));
    s = poly.toString();
    assert(s == "POLYGON ((0 1, 1 1, 1 0, 0 0, 0 1))", s); // same as polygon as hole is out of the polygon
}

@("Geometry - union")
@safe unittest
{
    auto p1 = Geometry.createPolygon(
        Point(0, 0),
        Point(1, 0),
        Point(1, 1),
        Point(0, 1),
        Point(0, 0)
    );

    auto p2 = Geometry.createPolygon(
        Point(2, 0),
        Point(3, 0),
        Point(3, 1),
        Point(2, 1),
        Point(2, 0)
    );

    auto p3 = Geometry.createPolygon(
        Point(1, 0.4),
        Point(2, 0.4),
        Point(2, 0.6),
        Point(1, 0.6),
        Point(1, 0.4)
    );

    {
        auto g = Geometry.createCollection(GEOSGeomTypes.GEOS_GEOMETRYCOLLECTION, p1.clone, p2.clone);
        assert(g.typeId == GEOSGeomTypes.GEOS_GEOMETRYCOLLECTION);

        auto u = g.union_();
        assert(g.typeId == GEOSGeomTypes.GEOS_GEOMETRYCOLLECTION); // nothing common between them
    }

    {
        auto g = Geometry.createCollection(GEOSGeomTypes.GEOS_GEOMETRYCOLLECTION, p1, p2, p3);
        assert(g.typeId == GEOSGeomTypes.GEOS_GEOMETRYCOLLECTION);

        auto u = g.union_();
        auto s = u.toString();
        assert(s == "POLYGON ((1 0, 0 0, 0 1, 1 1, 1 0.6, 2 0.6, 2 1, 3 1, 3 0, 2 0, 2 0.4, 1 0.4, 1 0))", s);
        assert(u.typeId == GEOSGeomTypes.GEOS_POLYGON); // merged to one polygon
    }
}

/**
 * Wrapper for `GEOSPreparedGeometry`.
 * It is created from normal `Geometry` and takes ownership of it (so use copy when original geometry is needed)
 */
struct PreparedGeometry
{
    private Geometry g;
    private const(GEOSPreparedGeometry)* pg;

    /// Constructor
    this(Geometry g)
    in (g.g, "Uninitialized geometry")
    {
        this.g = g.move();
        this.pg = GEOSPrepare_r(ctx, this.g.g);
        assert(pg !is null, "GEOSPrepare()");
    }

    ~this() @trusted nothrow @nogc
    {
        destroy(g);
        if (pg) GEOSPreparedGeom_destroy_r(ctx, pg);
    }

    /// Checks if the Geometry is initialized
    bool opCast(T)() const if (is(T == bool))
    {
        return pg !is null;
    }

    /// Access internally held geometry
    ref const(Geometry) geometry() return const @safe pure nothrow @nogc
    {
        return g;
    }
}

template checkRelationship(string rel)
{
    bool checkRelationship(GA, GB...)(auto ref const(GA) geomA, auto ref const(GB) geomB) @trusted nothrow @nogc
        if ((is(GA == Geometry) || is(GA == PreparedGeometry)) && GB.length >= 1)
    in (geomA, "Uninitialized geometry")
    {
        static if (is(GB[0] == Geometry)) {
            static assert(GB.length == 1, "Too many arguments");
            alias chg = geomB[0];
        }
        else static if (allSatisfy!(isNumeric, GB)) {
            static Geometry pt;
            alias chg = pt;
            static if (GB.length == 2) pt = Point(geomB[0], geomB[1]);
            else static if (GB.length == 3) pt = Point(geomB[0], geomB[1], geomB[2]);
            else static assert(0, "Too many coordinates for Point");
        }
        else static if (isPointStruct!(GB[0])) {
            static assert(GB.length == 1, "Too many arguments");
            static Geometry pt;
            alias chg = pt;
            pt = geomB[0];
        }
        else static assert(0, "Unsupported arguments combination");

        assert(chg, "Uninitialized geometry");
        static if (is(G == Geometry))
            auto r = mixin("GEOS" ~ rel ~ "_r")(ctx, geomA.g, chg.g);
        else
            auto r = mixin("GEOSPrepared" ~ rel ~ "_r")(ctx, geomA.pg, chg.g);
        assert(r != 2);
        return r == 1;
    }
}

/// Opposite to intersects
/// Overlaps, Touches, Within all imply geometries are not spatially disjoint. If any of the
/// aforementioned returns true, then the geometries are not spatially disjoint. Disjoint implies
/// false for spatial intersection.
alias disjoint = checkRelationship!"Disjoint";

/// Returns TRUE if A and B intersect, but their interiors do not intersect. Equivalently, A and B
/// have at least one point in common, and the common points lie in at least one boundary. For
/// Point/Point inputs the relationship is always FALSE, since points do not have a boundary.
alias touches = checkRelationship!"Touches";

/// Compares two geometries and returns true if they intersect. Geometries intersect if they have any point in common.
alias intersects = checkRelationship!"Intersects";

/// Compares two geometry objects and returns true if their intersection "spatially cross", that is,
/// the geometries have some, but not all interior points in common.
alias crosses = checkRelationship!"Crosses";

/// Tests if no points of A lie in the exterior of B, and A and B have at least one interior point in common.
alias within = checkRelationship!"Within";

/// Returns TRUE if geometry B is completely inside geometry A. A contains B if and only if no
/// points of B lie in the exterior of A, and at least one point of the interior of B lies in the
/// interior of A.
alias contains = checkRelationship!"Contains";

/// Returns TRUE if geometry A and B "spatially overlap". Two geometries overlap if they have the
/// same dimension, each has at least one point not shared by the other (or equivalently neither
/// covers the other), and the intersection of their interiors has the same dimension. The overlaps
/// relationship is symmetrical.
alias overlaps = checkRelationship!"Overlaps";

/// Tests the spatial equality of two geometries
alias equals = checkRelationship!"Equals";

/// Returns true if no point in Geometry/Geography B is outside Geometry/Geography A. Equivalently,
/// tests if every point of geometry B is inside (i.e. intersects the interior or boundary of)
/// geometry A.
alias covers = checkRelationship!"Covers";

/// Returns true if no point in Geometry/Geography A lies outside Geometry/Geography B.
/// Equivalently, tests if every point of geometry A is inside (i.e. intersects the interior or
/// boundary of) geometry B.
alias coveredBy = checkRelationship!"CoveredBy";

@("PreparedGeometry")
unittest
{
    auto poly = PreparedGeometry(Geometry.createPolygon(
        Point(0, 0),
        Point(1, 0),
        Point(1, 1),
        Point(0, 1),
        Point(0, 0)
    ));

    assert(poly.contains(0.5, 0.5));
    assert(poly.contains(Geometry.createPolygon(
        Point(0.4, 0.4),
        Point(0.6, 0.4),
        Point(0.6, 0.6),
        Point(0.4, 0.6),
        Point(0.4, 0.4)
    )));
}

/// Parameters to control buffered shape creation
struct BufferParameters
{
    GEOSBufCapStyles    endCapStyle = GEOSBufCapStyles.GEOSBUF_CAP_FLAT;
    GEOSBufJoinStyles   joinStyle   = GEOSBufJoinStyles.GEOSBUF_JOIN_ROUND;

    /// Set the mitre limit of a GEOSBufferParams to the desired size. For acute angles, a mitre
    /// join can extend very very far from the input geometry, which is probably not desired. The
    /// mitre limit places an upper bound on that.
    /// Note: Used only with mitre join style
    double mitreLimit = 5.0;

    /// Set the number of segments to use to stroke each quadrant of circular arcs generated by the
    /// buffering process. More segments means a smoother output, but with larger size.
    int quadSegments = 8;

    /// Sets whether the computed buffer should be single-sided. A single-sided buffer is
    /// constructed on only one side of each input line.
    ///
    /// The side used is determined by the sign of the buffer distance:
    ///   * a positive distance indicates the left-hand side
    ///   * a negative distance indicates the right-hand side
    bool singleSided;
}

/// Generates buffer from geometry.
/// Params:
///   - geom    = source geometry
///   - width   = the distance by which to expand the geometry (or contract if the value is negative) on each side.
///   - params  = parameters to control the shape of generated buffer
@trusted nothrow @nogc
Geometry buffer(ref const(Geometry) geom, double width, BufferParameters params = BufferParameters.init)
in (geom, "Uninitialized geometry")
{
    auto bp = GEOSBufferParams_create_r(ctx);
    if (!bp) onOutOfMemoryError();
    scope (exit) GEOSBufferParams_destroy_r(ctx, bp);

    auto r = GEOSBufferParams_setEndCapStyle_r(ctx, bp, params.endCapStyle);
    assert(r == 1, "GEOSBufferParams_setEndCapStyle()");

    r = GEOSBufferParams_setJoinStyle_r(ctx, bp, params.joinStyle);
    assert(r == 1, "GEOSBufferParams_setJoinStyle()");

    r = GEOSBufferParams_setMitreLimit_r(ctx, bp, params.mitreLimit);
    assert(r == 1, "GEOSBufferParams_setMitreLimit()");

    r = GEOSBufferParams_setQuadrantSegments_r(ctx, bp, params.quadSegments);
    assert(r == 1, "GEOSBufferParams_setQuadrantSegments()");

    r = GEOSBufferParams_setSingleSided_r(ctx, bp, params.singleSided);
    assert(r == 1, "GEOSBufferParams_setSingleSided()");

    Geometry ret;
    ret.g = GEOSBufferWithParams_r(ctx, geom.g, bp, width);
    if (!ret.g) assert(0, "GEOSBufferWithParams()");

    return ret;
}

@("Geometry buffer")
@safe unittest
{
    auto ln = Geometry.createLineString(Point(-1, 0.5), Point(1, 0.5));
    auto bln = ln.buffer(0.5);
    assert(bln.toString() == "POLYGON ((1 1, 1 0, -1 0, -1 1, 1 1))", bln.toString());
}

/// Corresponds with PostGIS output numbers
enum LineCrossDirection
{
    noCross = 0,                /// No crossing between lines
    left = -1,                  /// Line crosses the first one to the left
    right = 1,                  /// Line crosses the first one to the right
    multicrossLeft = -2,        /// Lines have multiple crosses, but at the end line is crossed to the left
    multicrossRight = 2,        /// Lines have multiple crosses, but at the end line is crossed to the right
    multicrossFirstLeft = -3,   /// Lines have multiple crosses with no side change, but first cross was to the left
    multicrossFirstRight = 3    /// Lines have multiple crosses with no side change, but first cross was to the right
}

/**
 * Determines direction in which the line A is crossed by line B.
 *
 * Only linestring geometries are valid as input.
 */
LineCrossDirection lineCrossingDirection(ref const Geometry lineA, ref const Geometry lineB) @safe nothrow @nogc
in (lineA.typeId == GEOSGeomTypes.GEOS_LINESTRING && lineB.typeId == GEOSGeomTypes.GEOS_LINESTRING, "Geometry must be of LINESTRING type")
{
    enum
    {
        SEG_NO_INTERSECTION = 0,
        SEG_COLINEAR = 1,
        SEG_CROSS_LEFT = 2,
        SEG_CROSS_RIGHT = 3
    }

    // check if segment envelopes has some overlapping
    static bool segInteract(const Point p1, const Point p2, const Point q1, const Point q2) nothrow @nogc @safe pure
    {
        pragma(inline, true);
        double minq = min(q1.x, q2.x);
        double maxq = max(q1.x, q2.x);
        double minp = min(p1.x, p2.x);
        double maxp = max(p1.x, p2.x);

        if (minp > maxq || maxp < minq)
            return false;

        minq = min(q1.y, q2.y);
        maxq = max(q1.y, q2.y);
        minp = min(p1.y, p2.y);
        maxp = max(p1.y, p2.y);

        if (minp > maxq || maxp < minq)
            return false;

        return true;
    }

    /*
     * Return -1  if point Q is left of segment P
     * Return  1  if point Q is right of segment P
     * Return  0  if point Q in on segment P
     */
    static int segmentSide(const Point p1, const Point p2, const Point q) @safe nothrow @nogc pure
    {
        pragma(inline, true);
        double side = ( (q.x - p1.x) * (p2.y - p1.y) - (p2.x - p1.x) * (q.y - p1.y) );
        return (side > 0) - (side < 0);
    }

    // returns the kind of crossing behavior of line segment 1 (constructed from p1 and p2) and line
    // segment 2 (constructed from q1 and q2)
    static int segmentIntersects(const Point p1, const Point p2, const Point q1, const Point q2) nothrow @nogc @safe pure
    {
        pragma(inline, true);

        // No envelope interaction => we are done
        if (!segInteract(p1, p2, q1, q2)) return SEG_NO_INTERSECTION;

        // Are the start and end points of q on the same side of p?
        auto pq1 = segmentSide(p1, p2, q1);
        auto pq2 = segmentSide(p1, p2, q2);
        if ((pq1>0 && pq2>0) || (pq1<0 && pq2<0)) return SEG_NO_INTERSECTION;

        // Are the start and end points of p on the same side of q?
        auto qp1 = segmentSide(q1, q2, p1);
        auto qp2 = segmentSide(q1, q2, p2);
        if ((qp1>0 && qp2>0) || (qp1<0 && qp2<0)) return SEG_NO_INTERSECTION;

        // Nobody is on one side or another? Must be colinear.
        if (pq1 == 0 && pq2 == 0 && qp1 == 0 && qp2 == 0) return SEG_COLINEAR;

        // Second point of p or q touches, it's not a crossing.
        if (pq2 == 0 || qp2 == 0) return SEG_NO_INTERSECTION;

        // First point of p touches, it's a "crossing"
        if (pq1 == 0)
            return pq2 > 0 ? SEG_CROSS_RIGHT : SEG_CROSS_LEFT;

        // First point of q touches, it's a crossing.
        if (qp1 == 0)
            return pq1 < pq2 ? SEG_CROSS_RIGHT : SEG_CROSS_LEFT;

        // The segments cross, what direction is the crossing?
        return pq1 < pq2 ? SEG_CROSS_RIGHT : SEG_CROSS_LEFT;
    }

    auto ptsA = lineA.coordSeq();
    auto ptsB = lineB.coordSeq();

    // One-point lines can't intersect (and shouldn't exist)
    if (ptsA.length < 2 || ptsB.length < 2) return LineCrossDirection.noCross;

	int leftCrosses = 0;
	int rightCrosses = 0;
	int firstCross = 0;

    Point p1, p2, q1, q2;
    q1 = ptsB[0];

    for (uint i = 1; i < ptsB.length; i++)
    {
        // Update second point of q to next value
        q2 = ptsB[i];

        /// Initialize first point of p
        p1 = ptsA[0];
        for (uint j = 1; j < ptsA.length; j++)
        {
            /// Update second point of p to next value
            p2 = ptsA[j];

            int cross = segmentIntersects(p1, p2, q1, q2);

            if (cross == SEG_CROSS_LEFT)
            {
                leftCrosses++;
                if (!firstCross)
                    firstCross = SEG_CROSS_LEFT;
            }
            else if (cross == SEG_CROSS_RIGHT)
            {
                rightCrosses++;
                if (!firstCross)
                    firstCross = SEG_CROSS_LEFT;
            }
            else if (cross == SEG_COLINEAR)
            {
                // TODO: Crossing at a co-linearity can be handled by extending
                // segment to next vertex and seeing if the end points straddle
                // the co-linear segment.
                // continue;
            }

            // Prepare new line segment start point
            p1 = p2;
        }

        // Prepare new line segment start point
        q1 = q2;
    }

    if (!leftCrosses && !rightCrosses) return LineCrossDirection.noCross;

    if (!leftCrosses  && rightCrosses == 1) return LineCrossDirection.right;
    if (!rightCrosses && leftCrosses  == 1) return LineCrossDirection.left;

    if (leftCrosses - rightCrosses > 0) return LineCrossDirection.multicrossLeft;
    if (leftCrosses - rightCrosses < 0) return LineCrossDirection.multicrossRight;

    if (leftCrosses - rightCrosses == 0 && firstCross == SEG_CROSS_LEFT) return LineCrossDirection.multicrossFirstLeft;
    if (leftCrosses - rightCrosses == 0 && firstCross == SEG_CROSS_RIGHT) return LineCrossDirection.multicrossFirstRight;

    return LineCrossDirection.noCross;
}

@("lineCrossingDirection - single segment")
@safe nothrow @nogc unittest
{
    // Vertical line segment from 0,0 to 0,1
    auto l1 = Geometry.createLineString(CoordSequence(Point(0,0), Point(0,1)));
    assert(lineCrossingDirection(l1, l1) == LineCrossDirection.noCross);

    // Horizontal line segment - crossing in the middle
    auto l2 = Geometry.createLineString(CoordSequence(Point(-0.5, 0.5), Point(0.5, 0.5)));
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.right);
    assert(lineCrossingDirection(l2, l1) == LineCrossDirection.left);

    // reverted direction
    l2 = Geometry.createLineString(CoordSequence(Point(0.5, 0.5), Point(-0.5, 0.5)));
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.left);
    assert(lineCrossingDirection(l2, l1) == LineCrossDirection.right);

    // Horizontal line segment - crossing at top end vertex (end crossings don't count)
    l2 = Geometry.createLineString(CoordSequence(Point(-0.5, 1), Point(0.5, 1)));
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.noCross);

    // Horizontal line segment - crossing at bottom end vertex
    l2 = Geometry.createLineString(CoordSequence(Point(-0.5, 0), Point(0.5, 0)));
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.right);

    // Horizontal line segment - no crossing
    l2 = Geometry.createLineString(CoordSequence(Point(-0.5, 2), Point(0.5, 2)));
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.noCross);

    // Vertical line segment - no crossing
    l2 = Geometry.createLineString(CoordSequence(Point(-0.5, 0), Point(-0.5, 1)));
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.noCross);
}

@("lineCrossingDirection - long lines")
@safe nothrow @nogc unittest
{
    // Vertical line with vertices at y integers
    auto l1 = Geometry("LINESTRING(0 0, 0 1, 0 2, 0 3, 0 4)");

    // Two crossings at segment midpoints
    auto l2 = Geometry("LINESTRING(1 1, -1 1.5, 1 3, 1 4, 1 5)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.multicrossFirstLeft);

    // One crossing at interior vertex
    l2 = Geometry("LINESTRING(1 1, 0 1, -1 1, -1 2, -1 3)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.left);

    // Two crossings at interior vertices
    l2 = Geometry("LINESTRING(1 1, 0 1, -1 1, 0 3, 1 3)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.multicrossFirstLeft);

    // Two crossings, one at the first vertex on at interior vertex
    l2 = Geometry("LINESTRING(1 0, 0 0, -1 1, 0 3, 1 3)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.multicrossFirstLeft);

    // Two crossings, one at the first vertex on the next interior vertex
    l2 = Geometry("LINESTRING(1 0, 0 0, -1 1, 0 1, 1 2)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.multicrossFirstLeft);

    // Three crossings, two at midpoints, one at vertex
    l2 = Geometry("LINESTRING(0.5 1, -1 0.5, 1 2, -1 2, -1 3)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.multicrossLeft);

    // One mid-point co-linear crossing
    l2 = Geometry("LINESTRING(1 1, 0 1.5, 0 2.5, -1 3, -1 4)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.left);

    // One on-vertices co-linear crossing
    l2 = Geometry("LINESTRING(1 1, 0 1, 0 2, -1 4, -1 4)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.left);

    // No crossing, but end on a co-linearity
    l2 = Geometry("LINESTRING(1 1, 1 2, 1 3, 0 3, 0 4)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.noCross);

    // Real life sample
    l1 = Geometry("LINESTRING(2.99 90.16,71 74,20 140,171 154)");
	l2 = Geometry("LINESTRING(25 169,89 114,40 70,86 43)");
    assert(lineCrossingDirection(l1, l2) == LineCrossDirection.multicrossRight);
}

private:
GEOSContextHandle_t ctx; // thread GEOS context
GEOSWKTWriter* wktWriter;
GEOSWKTReader* wktReader;
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

    // Workaround: above code doesn't work due to the https://issues.dlang.org/show_bug.cgi?id=21425
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
