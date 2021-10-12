/*
	��������� ��� ���������� ������� ������� ������ � ������� ������� ��� ���������� ����.
	��������� ������: shock_tube.dat
	������ ��������� ������:
		# ������������ ������� ����� 
		1.0	! ������ ������� �����
		0.5	! ������� �����������
		0.2	! ����� ������� �������
		2.0	! ��������� �����
		0.0	! �������� �����
		2.0	! �������� �����
		1.0	! ��������� ������
		0.0	! �������� ������
		1.0	! �������� ������
		1000	! ���������� �����
	�������� ������: out.dat (����������, ���������, ��������, ��������, ���������� �������)
	
	������ �������� � ������: 
	E.F. Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics". Springer-Verlag, 1997          
*/
#include <stdio.h>
#include <math.h>

#define EPS			1.0e-6	// �������� ������� ����������� ���������
#define MAX_ITER	20		// ���������� �������� ��� ������� ����������� ���������

/***** ��������� �� ���������� �������� *****/
double	dgamma = 1.4, 
		g1 = (dgamma - 1.0)/(2.0*dgamma),
		g2 = (dgamma + 1.0)/(2.0*dgamma),
		g3 = 2.0*dgamma/(dgamma - 1.0),
		g4 = 2.0/(dgamma - 1.0),
		g5 = 2.0/(dgamma + 1.0),
		g6 = (dgamma - 1.0)/(dgamma + 1.0),
		g7 = (dgamma - 1.0)/2.0,
		g8 = dgamma - 1.0;

/***** ��������� ������� ����� *****/
	// ��������� �����
double	dl,  // ���������
		ul,  // ��������
		pl,  // ��������
		cl;  // �������� �����

	// ��������� ������
double	dr,  // ���������
		ur,  // ��������
		pr,  // ��������
		cr;  // �������� �����

/**************************************************
 ��������� ���������, �������� � �������� � �����
**************************************************/
void sample(double &pm, double &um, double &s, 
			double &d, double &u, double &p)
{
	double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

	if(s <= um)
	{
		// ����� ����� �� ����������� �������
		if(pm <= pl)
		{
			// ����� ����� ����������
			shl = ul - cl;

			if(s <= shl)
			{
				// ����� �����
				d = dl;
				u = ul;
				p = pl;
			}
			else
			{
				cml = cl * pow(pm/pl,g1);
				stl = um - cml;

				if( s > stl )
				{
					// ����� ����� �� ����������� �������
					d = dl * pow(pm/pl,1.0/dgamma);
					u = um;
					p = pm;
				}
				else
				{
					// ����� � ����� ����������
					u = g5 * (cl + g7*ul + s);
					c = g5 * (cl + g7*(ul-s));
					d = dl*pow(c/cl,g4);
					p = pl*pow(c/cl,g3);
				}
			}
		}
		else
		{
			// ����� ������� �����
			pml = pm/pl;
			sl  = ul - cl * sqrt(g2*pml + g1);

			if( s <= sl )
			{
				// ����� ����� �� ������� �����
				d = dl;
				u = ul;
				p = pl;
			}
			else
			{
				// ����� ����� �� ����������� �������
				d = dl * (pml + g6)/(pml*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
	}
	else
	{
		// ����� ������ �� ����������� �������
		if( pm > pr )
		{
			// ������ ������� �����
			pmr = pm/pr;
			sr  = ur + cr * sqrt(g2*pmr + g1);

			if( s >= sr )
			{
				// ����� ������ �� ������� �����
				d = dr;
				u = ur;
				p = pr;
			}
			else
			{
				// ����� ������ �� ����������� �������
				d = dr * (pmr + g6)/(pmr*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
		else
		{
			// ������ ����� ����������
			shr = ur + cr;

			if( s >= shr )
			{
				// ����� ������ �� ����� ����������
				d = dr;
				u = ur;
				p = pr;
			}
			else
			{
				cmr = cr * pow(pm/pr,g1);
                str = um + cmr;

				if( s <= str )
				{
					// ����� ������ �� ����������� �������
					d = dr * pow(pm/pr,1.0/dgamma);
					u = um;
					p = pm;
				}
				else
				{
					// ����� � ����� ����� ����������
					u = g5 * (-cr + g7*ur + s);
					c = g5 * (cr - g7*(ur-s));
					d = dr * pow(c/cr,g4);
					p = pr * pow(c/cr,g3);
				}
			}
		}
	}

}

/**************************************************
		������� ��� ����� �� ������
**************************************************/
void prefun(double &f, double &fd, double &p, 
			double &dk, double &pk, double &ck)
{
	double ak, bk, pratio, qrt; 
	
	if(p <= pk)
	{
		// ����� ����������
		pratio = p/pk;
		f	   = g4*ck*(pow(pratio,g1) - 1.0);
		fd     = (1.0/(dk*ck))*pow(pratio,-g2);

	}
	else
	{
		// ������� �����
		ak  = g5/dk;
		bk  = g6*pk;
		qrt = sqrt(ak/(bk + p));
		f   = (p - pk)*qrt;
		fd  = (1.0 - 0.5*(p - pk)/(bk + p)) * qrt;
	}

}

/**************************************************
		���������� ���������� �����������
**************************************************/
double guessp()
{
	double	cup, gel, ger, 
			pmax, pmin, ppv, pq, pm, 
			ptl, ptr, 
			qmax, quser, um;

	/*** ���������� ����������� �������� �� PVRS ������� ������ ***/
	quser = 2.0;
	cup = 0.25 * (dl + dr) * (cl + cr);
	ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
	ppv = ppv > 0.0 ? ppv : 0.0;
	pmin = pl > pr ? pr : pl;
	pmax = pl > pr ? pl : pr;
	qmax = pmax/pmin;

    if( qmax<=quser && pmin<=ppv && ppv<=pmax )
	{
		pm = ppv;
	}
	else
	{
		if(ppv < pmin)
		{
			// ��� ����� ����������
            pq  = pow(pl/pr,g1);
            um  = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
            ptl = 1.0 + g7*(ul - um)/cl;
            ptr = 1.0 + g7*(um - ur)/cr;
            pm  = 0.5*(pl*pow(ptl,g3) + pr*pow(ptr,g3));
		}
		else
		{
			// ��� ������� �����
			gel = sqrt((g5/dl)/(g6*pl + ppv));
			ger = sqrt((g5/dr)/(g6*pr + ppv));
            pm  = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
		}

	}
	return pm;
}

/************************************************** 
����������� �������� � �������� ����������� �������			 
**************************************************/
void starpu(double &p, double &u)
{
	int i;
	
	double pstart,	// ��������� ����������� �������� 
		   pold,	// ���������� ����������� ��������
		   udiff,	// �������� ���������
		   change;	// �������� ��������

	double	fl, fld, // ������� ��������
			fr, frd;	

	/*** ���������� ���������� ����������� ***/
	pstart = guessp();

	/*** ���������� ����������� ***/
	pold = pstart;

	/*** �������� ��������� ***/
	udiff = ur - ul;

	/*** ����� ������� ��� ����������� �������� ***/

	// �������� ����� ������� ������������� ��������
	change = 10.0 * EPS;

	for(i=0 ; i<MAX_ITER && change>EPS ; i++)
	{
		// ������� ��� ����� �����
		prefun(fl,fld,pold,dl,pl,cl);

		// ������� ��� ������ �����
		prefun(fr,frd,pold,dr,pr,cr);

		// ��������� ����������� ��������
		p = pold - (fl + fr + udiff) / (fld + frd);

		// �������� ����� ������� ������������� ��������
		change = 2.0 * fabs((p-pold)/(p+pold));

		// ���� �������� �������������, �� �������� ���
		if(p < 0.0) p =0.0;

		pold = p;
	}

	// ����������� ��������
	u = 0.5 * (ul + ur + fr - fl);
}

/************************************************** 
				�������� ������� 
**************************************************/
int main()
{
	char string[128];	// ������ � ����� ������������ ������� �����

	double	diaph,	// ������� �����������
			len,	// ������ ������� [0;len]
			time,	// ����� ������� �������
			pm,		// �������� �� �������� ����������� �������
			um,		// �������� ����������� �������
			dx,		// ��� �� ������������
			xpos,	// ������� ����� �������
			ds,		// ��������� � ����� �������
			ps,		// �������� � ����� �������
			us,		// �������� � ����� �������
			s;		// ��������������

	int numcells;	// ���������� �����
	int i;			// ����� ������

	FILE *fout,		// ���� �����������
		 *fin;		// ���� � �������

	//*** ������� ���������� ����� ***
		fin = fopen("shock_tube.dat","r");
		fgets(string,sizeof(string),fin);
	// ������ ����� [0;len]
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&len);
	// ������� �����������
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&diaph);
	// ����� ������� �������
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&time);
	// ��������� �����
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&dl);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&ul);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&pl);

		cl	= sqrt(dgamma*pl/dl);
	// ��������� ������
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&dr);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&ur);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&pr);

		cr	= sqrt(dgamma*pr/dr);
	// ���������� �����
		fgets(string,sizeof(string),fin);
		sscanf(string,"%d",&numcells);

	//********************************


	/*** ����������� ������� ������� ***/
	if( g4*(cl+cr) <= (ur-ul) )
	{
		printf("*** Generated vacuum ***\n");
		return -1;
	}

	/*** ���������� ������� ������� �������� � �������� ����������� ������� ***/
	starpu(pm,um); 

	/*** ����������� ���� (������� ������) ***/
	dx = len/double(numcells);

	/*** ����� ������� ***/
	fout = fopen("out.dat","w");
	for(i=1 ; i<=numcells ; i++)
	{
		xpos = i*dx - 0.5*dx;
		s = (xpos - diaph)/time;
		sample(pm,um,s,ds,us,ps);
		/***************************************************** 
		 | x | density | speed | pressure | internal energy |
		*****************************************************/
		//fprintf(fout,"%lf %lf %lf %lf %lf\n",xpos,ds,us,ps,ps/ds/g8);
		fprintf(fout,"%lf %lf %lf %lf\n",xpos,ds,us,ps);
	}
	fclose(fout);
	return 0;
}
