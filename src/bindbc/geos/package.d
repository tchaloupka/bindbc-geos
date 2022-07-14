module bindbc.geos;

public import bindbc.geos.wrapper;

version (unittest):

import bindbc.geos.libgeos;
import core.stdc.stdarg;
import core.stdc.stdio;

extern (C) nothrow @nogc
void msgHandler(const(char)* fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
}

// Basic library usage test from [building-a-program](https://libgeos.org/usage/c_api/#building-a-program) example.
unittest
{
    printf("Using GEOS Version: %s\n", GEOSversion());
    initGEOS(&msgHandler, &msgHandler);
    scope (exit) bindbc.geos.libgeos.finishGEOS();

    /* Read WKT into geometry object */
    GEOSWKTReader* reader = GEOSWKTReader_create();
    GEOSGeometry* geom_a = GEOSWKTReader_read(reader, "POINT(12 34)");

    /* Convert result to WKT */
    GEOSWKTWriter* writer = GEOSWKTWriter_create();
    char* wkt = GEOSWKTWriter_write(writer, geom_a);
    printf("Geometry WKT: %s\n", wkt);

    static if (geosSupport >= GEOSSupport.geos_3_10) {
        /* Convert result to GeoJSON */
        GEOSGeoJSONWriter* jwriter = GEOSGeoJSONWriter_create();
        wkt = GEOSGeoJSONWriter_writeGeometry(jwriter, geom_a, -1);
        printf("Geometry GeoJSON: %s\n", wkt);
        GEOSGeoJSONWriter_destroy(jwriter);
    }

    /* Clean up allocated objects */
    GEOSWKTReader_destroy(reader);
    GEOSWKTWriter_destroy(writer);
    GEOSGeom_destroy(geom_a);
    GEOSFree(wkt);
}
