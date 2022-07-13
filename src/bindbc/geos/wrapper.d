/**
 * More D idiomatic API wrapper around GEOS C API
 */
module bindbc.geos.wrapper;

import bindbc.geos.libgeos;

@safe nothrow @nogc:

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
    GEOSGeomTypes typeId() const nothrow @nogc @trusted
    in (g !is null, "Can't get type of uninitialized geometry")
    {
        auto tid = GEOSGeomTypeId(g);
        assert(tid >= 0, "Error getting geometry typeId");
        return cast(GEOSGeomTypes)tid;
    }

    /**
     * Create a new copy of the geometry.
     */
    Geometry clone() const nothrow @nogc @trusted
    {
        Geometry c;
        if (this.g) c.g = GEOSGeom_clone(this.g);
        return c;
    }

    @disable this(this); /// Copy is not allowed - use move() or clone()

    ~this() nothrow @nogc @trusted
    {
        if (g) GEOSGeom_destroy(g);
    }

    private:
    GEOSGeometry* g;
}
