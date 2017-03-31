#include "medit.h"
#include "extern.h"
#include "sproto.h"

extern int getIso(pScene sc,double norm,double *hsv);


GLuint geomList(pScene sc,pMesh mesh) {
  GLuint     list = 0;
  pMaterial  pm;
  pEdge      pr;
  pPoint     ppt,pp0,pp1;
	pSolution  ps0;
  double     dd,rgb[3],val;
  float      n[3];
  int        k,it = 0,nm;
	char       nodata;
  static float green[4] = {0.0, 1.0, 0.0, 1.0};
  static float rouge[4] = {1.0, 0.0, 0.0, 1.0};
  static float jaune[4] = {1.0, 1.0, 0.0, 1.0};
  static float orange[4]= {1.0, 0.65, 0.0, 1.0};
  static double hsv[3] = { 0.0, 1.0, 0.8 };
  
  /* default */
  if ( mesh->na+mesh->nc+mesh->np == 0 )  return(0);
  nodata = egal(sc->iso.val[0],sc->iso.val[MAXISO-1]);

  /* create display list */
  list = glGenLists(1);
  if ( !list )  return(0);
  glNewList(list,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* draw corners, ridges and required items */
  if ( ddebug ) printf("construct point list\n");
  if ( mesh->ne ) {
    /*--- ancienne partie --- */
    glPointSize(sc->par.pointsize);
    glBegin(GL_POINTS);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( ppt->tag & M_UNUSED )
        continue;
      else if ( ppt->tag & M_CORNER )
        glColor3fv(rouge);
      else if ( ppt->tag & M_REQUIRED )
        glColor3fv(green);
      else if ( !(ppt->tag & M_RIDGE) )
        continue;
      else if ( sc->par.linc == 1 )
        glColor3fv(sc->par.edge);
      else
        glColor3fv(orange);
      it++;
      glVertex3f(ppt->c[0],ppt->c[1],ppt->c[2]);
    }
    glEnd();
    /*-------------------------*/
    glPointSize(1);
  }
  else {
		if ( nodata )
			glPointSize(sc->par.pointsize);
		else
			glPointSize(2.0);
    glBegin(GL_POINTS);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      n[0] = ppt->c[0] - sc->cx;
      n[1] = ppt->c[1] - sc->cy;
      n[2] = ppt->c[2] - sc->cz;
      dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd > 0.0 ) {
        dd = 1.0 / sqrt(dd);
        n[0] *= dd;
        n[1] *= dd;
        n[2] *= dd;
      }
	    if ( !nodata ) {
	      ps0 = &mesh->sol[k];
				val = ps0->bb;
				getIso(sc,val,hsv);
		    hsvrgb(hsv,rgb);
		    glColor3dv(rgb);
	      glVertex3f(ppt->c[0],ppt->c[1],ppt->c[2]);
			}
	    else {
				glColor3fv(rouge);
        glNormal3fv(n);
	      glVertex3f(ppt->c[0],ppt->c[1],ppt->c[2]);
			}
    }
		glEnd();

    /*if ( !nodata) {
		  glPointSize(sc->par.pointsize);
      glBegin(GL_POINTS);
		  for (k=1; k<=mesh->np; k++) {
			  ppt = &mesh->point[k];
			  if ( !ppt->ref )  continue;
	      ps0 = &mesh->sol[k];
				val = ps0->bb;
				getIso(sc,val,hsv);
		    hsvrgb(hsv,rgb);
		    glColor3dv(rgb);
			  glVertex3f(ppt->c[0],ppt->c[1],ppt->c[2]);
		  }
	    glEnd();
		}*/
    it = mesh->np;
  }

  /* draw edges */
  if ( ddebug )  printf("construct edge list\n");
  glLineWidth(sc->par.linewidth);
  glBegin(GL_LINES);
  for (k=1; k<=mesh->na; k++) {
    pr = &mesh->edge[k];
    if ( pr->v[0] > mesh->np || pr->v[1] > mesh->np )  continue;

    if ( pr->tag & M_RIDGE ) {
      if ( pr->tag & M_TAG )
	      glColor3fv(jaune);  /* ridge + ref en jaune */
      else
	      glColor3fv(rouge);  /* ridges en rouge */
    }
    else if ( !pr->ref ) {
      glColor3fv(sc->par.edge);
    }
    else {
      nm = matRef(sc,pr->ref);
      pm = &sc->material[nm];
      glColor3fv(pm->dif);
    }
    if ( sc->par.linc == 1 )  glColor3fv(sc->par.edge);
    pp0 = &mesh->point[pr->v[0]];
    pp1 = &mesh->point[pr->v[1]];
    glVertex3f(pp0->c[0],pp0->c[1],pp0->c[2]);
    glVertex3f(pp1->c[0],pp1->c[1],pp1->c[2]);
    it++;
  }
  glEnd();
  glLineWidth(1.0);
  glEndList();

  if ( it == 0 ) {
    glDeleteLists(list,1);
    return(0);
  }
  else
    return(list);
}

