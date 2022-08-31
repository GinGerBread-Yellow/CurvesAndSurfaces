#ifdef __APPLE__
#include <OpenGL/gl.h>
/* Just in case we need these later
// References:
// http://alumni.cs.ucsb.edu/~wombatty/tutorials/opengl_mac_osx.html
// # include <OpenGL/gl.h>
// # include <OpenGL/glu.h>
*/
#else
#include <GL/gl.h>
#endif

#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

    
}
    

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }

    // TODO:
    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Make sure that this function computes all the appropriate
    // Vector3fs for each CurvePoint: V,T,N,B.
    // [NBT] should be unit and orthogonal.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.

    // cerr << "\t>>> evalBezier has been called with the following input:" << endl;

    // cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    // for( unsigned i = 0; i < P.size(); ++i )
    // {
    //     cerr << "\t>>> " << P[i] << endl;
    // }

    // cerr << "\t>>> Steps (type steps): " << steps << endl;
    // cerr << "\t>>> Returning empty curve." << endl;

    int pieces = P.size() / 3, end = P.size()-1;
    Curve ret(pieces*steps+1);

    
    // other way is to use BGT

    // B
    Matrix4f bernstein(1.0f, -3.0f, 3.0f, -1.0f,
        0.0f, 3.0f, -6.0f, 3.0f,
        0.0f, 0.0f, 3.0f, -3.0f,
        0.0f, 0.0f, 0.0f, 1.0f);
    // derivative
    // Matrix4f prime(0.0f, -3.0f, 6.0f, -3.0f,
    //     0.0f, 3.0f, -12.0f, 9.0f,
    //     0.0f, 0.0f, 6.0f, -9.0f,
    //     0.0f, 0.0f, 0.0f, 3.0f);
    Matrix4f prime(-3.f, 6.f, -3.f, 0.f, 
        3.f, -12.f, 9.f, 0.f, 
        0.f, 6.f, -9.f , 0.f, 
        0.f, 0.f, 3.f, 0.f);

    Vector3f bnorm(0.0f, 0.0f, 1.0f);
    int idx = 0;
    for(int j=0; j<end; j+=3) {
        Matrix4f control(P[j][0], P[j+1][0], P[j+2][0], P[j+3][0],
            P[j][1], P[j+1][1], P[j+2][1], P[j+3][1],
            P[j][2], P[j+1][2], P[j+2][2], P[j+3][2],
            0.0f, 0.0f, 0.0f, 0.0f);

        for(int i=0; i<steps; i++) {
            float t = (float)i/steps;

            Vector4f monomial(1, t, t*t, t*t*t);
            ret[idx].V = (control*bernstein*monomial).xyz();
            ret[idx].T = (control*prime*monomial).xyz().normalized();
            ret[idx].N = Vector3f::cross(bnorm, ret[idx].T).normalized();
            ret[idx].B = Vector3f::cross(ret[idx].T, ret[idx].N).normalized();
            bnorm = ret[idx].B;
            idx++;
        }
    }

    /**
    // Monomial basis vector
    Vector3f bnorm(0.0f, 0.0f, 1.0f);

    int idx = 0;
    for(int j=0; j<end; j+=3) {
        for(int i=0; i<steps; i++) {
            float alpha = (float)i/steps;
            Vector3f vp1 = Vector3f::lerp(P[j],P[j+1], alpha);
            Vector3f vp2 = Vector3f::lerp(P[j+1], P[j+2], alpha);
            Vector3f vp3 = Vector3f::lerp(P[j+2], P[j+3], alpha);
            Vector3f vp4 = Vector3f::lerp(vp1, vp2, alpha);
            Vector3f vp5 = Vector3f::lerp(vp2, vp3, alpha);
            Vector3f vp6 = Vector3f::lerp(vp4, vp5, alpha);
            ret[idx].V = vp6;
            ret[idx].T = (vp5-vp6).normalized();
            ret[idx].N = Vector3f::cross(bnorm, ret[idx].T).normalized();
            ret[idx].B = Vector3f::cross(ret[idx].T, ret[idx].N).normalized();
            bnorm = ret[idx].B;
            idx++;
            
        }
    }
    */

    // last node
    int last = pieces*steps;
    ret[last].V = P[end];
    ret[last].T = (P[end]-P[end-1]).normalized();
    ret[last].N = Vector3f::cross(bnorm, ret[last].T).normalized();
    ret[last].B = Vector3f::cross(ret[last].T, ret[last].N);
    // Right now this will just return this empty curve.
    return ret;

}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.

    // cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    // cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    // for( unsigned i = 0; i < P.size(); ++i )
    // {
    //     cerr << "\t>>> " << P[i] << endl;
    // }

    // cerr << "\t>>> Steps (type steps): " << steps << endl;
    // cerr << "\t>>> Returning empty curve." << endl;

    // Return an empty curve right now.

    int nnodes = 3 * P.size() - 8;
    vector< Vector3f > bezier(nnodes);
    int i=0;
    for(int seg = 1; seg < P.size()-1; seg++) {
        // 3 point
        Vector3f st = Vector3f::lerp(P[seg-1], P[seg], (float)2/3);
        Vector3f end = Vector3f::lerp(P[seg], P[seg+1], (float)1/3);
        Vector3f mid = Vector3f::lerp(st, end, (float)1/2);

        if(i>0) bezier[i++] = st;
        bezier[i++] = mid;
        if(i<nnodes) bezier[i++] = end;
    }
    return evalBezier(bezier, steps);
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

