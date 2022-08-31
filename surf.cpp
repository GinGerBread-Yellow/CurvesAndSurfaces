#include "surf.h"
#include "extra.h"
using namespace std;

namespace
{
    
    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i=0; i<profile.size(); i++)
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0){
                cerr<<"index "<<i
                    << ", V.z="<< profile[i].V[2] 
                    << ", T.z=" << profile[i].T[2] 
                    << ", N.z=" << profile[i].N[2] << '\n';
                return false;
        }
    
        return true;
    }
}

void meshing(Surface *surface, int nrow, int ncol) {
    vector< Tup3u > &vf = surface->VF;

    int k = ncol*nrow-ncol;
    vf.resize(2*(nrow-1)*(ncol-1));

    int idx = 0;
    for(unsigned i = 0; i < k; i++) {
        if((i+1)%ncol == 0) continue;

        vf[idx][0] = i+1;
        vf[idx][1] = i;
        vf[idx][2] = i+ncol;
        idx++;
        vf[idx][0] = i+1;
        vf[idx][1] = i+ncol;
        vf[idx][2] = i+ncol+1;
        idx++;
        
    }
    return;
}

Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    Surface surface;
    
    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.

    int n = profile.size();

    vector< Vector3f > &vv = surface.VV;
    vector< Vector3f > &vn = surface.VN;
    vector< Tup3u > &vf = surface.VF;
    vv.resize(n*(steps+1));
    vn.resize(n*(steps+1));


    unsigned idx = 0;
    for( unsigned v = 0; v < n; ++v) {

        for( unsigned i = 0; i <= steps; ++i )
        {
            // step from 0 to 2pi
            float t = 2.0f * M_PI * float( i ) / steps;

            Matrix3f rot = Matrix3f::rotateY(t);
            vv[idx] = rot * profile[v].V;
            vn[idx] = rot * profile[v].N;
            vn[idx].negate();
            idx++;
        }
    }

    meshing(&surface, n, steps+1);
    // idx = 0;
    // vf.resize((n-1)*(steps)*2);
    // for(unsigned v = 0; v < n-1; ++v ) {

    //     for( unsigned i = 0; i < steps; ++i ) {
    //         // find polygon vv[v*(steps+1)+i]
    //         unsigned l1 = v*(steps+1)+i, r1 = l1+1;
    //         unsigned l2 = (v+1)*(steps+1)+i, r2 = l2+1;

    //         vf[idx][0] = l1;
    //         vf[idx][1] = l2;
    //         vf[idx][2] = r1;
    //         idx++;
    //         vf[idx][0] = l2;
    //         vf[idx][1] = r2;
    //         vf[idx][2] = r1;
    //         idx++;
    //     }
    // }
    // cerr << "\t>>> makeSurfRev called (but not implemented).\n\t>>> Returning empty surface." << endl;
 
    return surface;
}
   
Surface makeGenCyl(const Curve &profile, const Curve &sweep )
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.
    int m = sweep.size(), n = profile.size();
    vector< Vector3f > &vv = surface.VV;
    vector< Vector3f > &vn = surface.VN;
    vector< Tup3u > &vf = surface.VF;
    vv.resize(m*n);
    vn.resize(m*n);

    int idx = 0;
    for(int j=0; j<n; ++j) {

        Vector4f pos(profile[j].V, 1.f);
        Vector4f norm(profile[j].N, 0.f);

        for(int i=0; i<m; ++i){
            const Vector3f &N = sweep[i].N;
            const Vector3f &B = sweep[i].B;
            const Vector3f &T = sweep[i].T;
            const Vector3f &V = sweep[i].V;
            Matrix4f transform(N[0], B[0], T[0], V[0],
                N[1], B[1], T[1], V[1],
                N[2], B[2], T[2], V[2],
                0.f, 0.f, 0.f, 1.f);

            vv[idx] = (transform*pos).xyz();
            vn[idx] = (transform*norm).xyz().normalized();
            vn[idx].negate();
            idx++;

        }
    }

    meshing(&surface, n, m);
    // cerr << "\t>>> makeGenCyl called (but not implemented).\n\t>>> Returning empty surface." <<endl;

    return surface;
}

void drawSurface(const Surface &surface, bool shaded)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    if (shaded)
    {
        // This will use the current material color and light
        // positions.  Just set these in drawScene();
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // This tells openGL to *not* draw backwards-facing triangles.
        // This is more efficient, and in addition it will help you
        // make sure that your triangles are drawn in the right order.
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    else
    {        
        glDisable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        glColor4f(0.4f,0.4f,0.4f,1.f);
        glLineWidth(1);
    }

    glBegin(GL_TRIANGLES);
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        glNormal(surface.VN[surface.VF[i][0]]);
        glVertex(surface.VV[surface.VF[i][0]]);
        glNormal(surface.VN[surface.VF[i][1]]);
        glVertex(surface.VV[surface.VF[i][1]]);
        glNormal(surface.VN[surface.VF[i][2]]);
        glVertex(surface.VV[surface.VF[i][2]]);
    }
    glEnd();

    glPopAttrib();
}

void drawNormals(const Surface &surface, float len)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_LIGHTING);
    glColor4f(0,1,1,1);
    glLineWidth(1);

    glBegin(GL_LINES);
    for (unsigned i=0; i<surface.VV.size(); i++)
    {
        glVertex(surface.VV[i]);
        glVertex(surface.VV[i] + surface.VN[i] * len);
    }
    glEnd();

    glPopAttrib();
}

void outputObjFile(ostream &out, const Surface &surface)
{
    
    for (unsigned i=0; i<surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (unsigned i=0; i<surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;
    
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j=0; j<3; j++)
        {
            unsigned a = surface.VF[i][j]+1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
