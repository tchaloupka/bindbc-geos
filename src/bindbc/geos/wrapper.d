/**
 * More D idiomatic API wrapper around GEOS C API
 */
module bindbc.geos.wrapper;

import core.stdc.stdarg;
import core.stdc.stdio;
import core.stdc.stdlib;
import bindbc.geos.libgeos;

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
    GEOSContext_setNoticeHandler_r(ctx, &msgHandlerNotice);
    GEOSContext_setErrorHandler_r(ctx, &msgHandlerError);
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

/**
 * Wrapper around GEOS library `GEOSGeometry*` type.
 * Geometry instance from GEOS is automatically freed using RAII.
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

    alias xMin = getPointCoordinate!GEOSGeom_getXMin_r; /// Finds the minimum X value in the geometry.
    alias xMax = getPointCoordinate!GEOSGeom_getXMax_r; /// Finds the maximum X value in the geometry.
    alias yMin = getPointCoordinate!GEOSGeom_getYMin_r; /// Finds the minimum Y value in the geometry.
    alias yMax = getPointCoordinate!GEOSGeom_getYMax_r; /// Finds the maximum Y value in the geometry.

    /// Get the total number of points in a geometry, of any type.
    int length() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    {
        return GEOSGetNumCoordinates_r(ctx, g);
    }

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
            xmin = this.xMin();
            ymin = this.yMin();
            xmax = this.xMax();
            ymax = this.yMax();
        }
    }

    /// ditto
    Extent extent() const @trusted nothrow @nogc
    {
        Extent ext;
        extent(ext.xmin, ext.ymin, ext.xmax, ext.ymax);
        return ext;
    }

    @disable this(this); /// Copy is not allowed - use move() or clone()

    ~this() @trusted nothrow @nogc
    {
        if (g) GEOSGeom_destroy_r(ctx, g);
    }

    /// Creates an empty point.
    static Geometry emptyPoint() @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    {
        Geometry p;
        p.g = GEOSGeom_createEmptyPoint_r(ctx);
        if (p.g is null) assert(0, "Failed to create emptyPoint");
        return p;
    }

    /// Creates a point geometry from a pair of coordinates.
    static Geometry createPoint(double x, double y) @trusted nothrow @nogc
    in (ctx !is null, "GEOS thread context not initialized")
    {
        Geometry p;
        p.g = GEOSGeom_createPointFromXY_r(ctx, x, y);
        if (p.g is null) assert(0, "Failed to create point");
        return p;
    }

    /// Generates geometry as [WKT](https://libgeos.org/specifications/wkt/) string
    void toString(S)(auto ref S sink) const @trusted
    {
        import core.stdc.string : strlen;
        char* s = GEOSWKTWriter_write_r(ctx, wktWriter, g);
        if (s is null) assert(0, "Failed to generate WKT from geometry");
        sink(s[0..strlen(s)]);
        GEOSFree_r(ctx, s);
    }

    /// ditto
    string toString() const @trusted
    {
        import core.stdc.string : strlen;
        char* s = GEOSWKTWriter_write_r(ctx, wktWriter, g);
        if (s is null) assert(0, "Failed to generate WKT from geometry");
        string str = s[0..strlen(s)].idup;
        GEOSFree_r(ctx, s);
        return str;
    }

    private:
    GEOSGeometry* g;

    double getPointCoordinate(alias fn)() const @trusted nothrow @nogc
    in (g !is null, "unitialized geometry")
    in (typeId() == GEOSGeomTypes.GEOS_POINT, "Operation is valid only on Point geometry")
    {
        double coord;
        auto r = fn(ctx, g, &coord);
        assert(r == 1, "Error getting coordinate");
        return coord;
    }
}

struct Extent
{
    double xmin;
    double ymin;
    double xmax;
    double ymax;
}

version (unittest)
{
    static this() @safe nothrow @nogc { initGEOS(); }
    static ~this() @safe nothrow @nogc { finishGEOS(); }
}

@safe
unittest
{
    import std.array : Appender;
    import std.math : isNaN;

    auto point = Geometry.createPoint(1, 2);
    assert(point.g !is null);
    assert(point.typeId == GEOSGeomTypes.GEOS_POINT);
    auto s = point.toString();
    assert(s == "POINT (1 2)", s);

    assert(point.x == 1);
    assert(point.y == 2);
    assert(point.z.isNaN);
    assert(point.length == 1);
    assert(point.dimensions == 0);
    assert(point.xMin == 1);
    assert(point.xMax == 1);
    assert(point.yMin == 2);
    assert(point.yMax == 2);

    auto ext = point.extent();
    assert(ext.xmin == 1);
    assert(ext.xmax == 1);
    assert(ext.ymin == 2);
    assert(ext.ymax == 2);

    Appender!string buf;
    point.clone().toString((const(char)[] s) {buf ~= s; });
    assert(buf.data == "POINT (1 2)", buf.data);
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
void msgHandlerNotice(const(char)* fmt, ...)
{

    va_list ap;
    va_start(ap, fmt);
    auto msg = makeMessage(fmt, ap);
    msgHandler(MessageLevel.notice, msg);
    free(msg.ptr);
    va_end(ap);
}

extern (C) nothrow @trusted
void msgHandlerError(const(char)* fmt, ...)
{

    va_list ap;
    va_start(ap, fmt);
    auto msg = makeMessage(fmt, ap);
    msgHandler(MessageLevel.error, msg);
    free(msg.ptr);
    va_end(ap);
}

// assembles the formatted message (needs to be freed afterwards)
char[] makeMessage(const char *fmt, va_list ap) @trusted nothrow @nogc
{
    /* Determine required size. */
    auto n = vsnprintf(null, 0, fmt, ap);

    if (n < 0)
        return null;

    size_t size = n + 1;      /* One extra byte for '\0' */
    auto p = cast(char*)malloc(size);
    if (p is null) return null;

    n = vsnprintf(p, size, fmt, ap);

    if (n < 0) {
        free(p);
        return null;
    }

    return p[0..n];
}
