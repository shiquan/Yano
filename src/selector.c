#include <R.h>
#include <Rinternals.h>
#include <R_ext/GraphicsEngine.h>

// https://stackoverflow.com/questions/217578
// https://wrfranklin.org/Research/Short_Notes/pnpoly.html
Rboolean pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
    int i, j;
    Rboolean c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
             (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
            c = !c;
  }
  return c;
}
SEXP points_in_polygen(SEXP x1, SEXP y1, SEXP x2, SEXP y2)
{
    if (length(x2) != length(y2) || length(x1) != length(y1)) error("unequal length");

    int i;
    int l = length(x2);
    SEXP ret;
    PROTECT(ret = allocVector(LGLSXP, l));

    int m = length(x1);
    double *x = malloc(m *sizeof(double));
    double *y = malloc(m *sizeof(double));

    for (i = 0; i < m; ++i) {
        x[i] = REAL(x1)[i];
        y[i] = REAL(y1)[i];
    }
    
    for (i = 0; i < l; ++i) {
        LOGICAL(ret)[i] = pnpoly(m, x, y, REAL(x2)[i], REAL(y2)[i]);
    }

    free(x);
    free(y);
    UNPROTECT(1);
    return ret;
}
SEXP selector()
{
    int n = 0;
    int m = 4;
    double *x = malloc(m *sizeof(double));
    double *y = malloc(m *sizeof(double));
    
    R_GE_gcontext gc;
    gc.lwd = 1;
    gc.lty = LTY_SOLID;
    gc.col = R_RGB(0,0,255);
    
    GEDevDesc *dd = GEcurrentDevice();

    if (!dd) {
        error("No active graphics device.");
    }

    if (!dd->dev->locator) {
            error("No locator complied.");
    }

    for (;;) {

        if (n == m) {
            m = m*2;
            x = realloc(x, m*sizeof(double));
            y = realloc(y, m*sizeof(double));
        }

    
        dd->dev->mode(1, dd->dev);
        if (!dd->dev->locator(x+n, y+n, dd->dev)) break;
    
        if (n > 0) {
            dd->dev->line(x[n-1], y[n-1], x[n], y[n], &gc, dd->dev);
        }
        dd->dev->mode(0, dd->dev);
        // if (dd->dev->circle) dd->dev->circle(x[n], y[n], 1.0, &gc, dd->dev);
        
        n++;
    }

    if (n > 2) {
        dd->dev->mode(1, dd->dev);
        dd->dev->line(x[n-1], y[n-1], x[0], y[0], &gc, dd->dev);
    }
    
    dd->dev->mode(0, dd->dev);

    SEXP ans, x0, y0, n0;
    PROTECT(ans = allocList(3));
    PROTECT(x0 = allocVector(REALSXP, n));
    PROTECT(y0 = allocVector(REALSXP, n));
    PROTECT(n0 = allocVector(INTSXP, 1));

    int i;
    for (i = 0; i < n; ++i) {
        REAL(x0)[i] = fromDeviceX(x[i], GE_INCHES, dd);
        REAL(y0)[i] = fromDeviceY(y[i], GE_INCHES, dd);
    }
    free(x);
    free(y);
    INTEGER(n0)[0] = n;
    SETCAR(ans, x0);
    SETCADR(ans, y0);
    SETCADDR(ans, n0);
    UNPROTECT(4);
    return ans;
}
