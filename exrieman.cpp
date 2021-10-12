/*
	Программа для нахождения точного решения задачи о распаде разрыва для идеального газа.
	Начальные данные: shock_tube.dat
	Пример начальных данных:
		# Конфигурация ударной трубы 
		1.0	! размер ударной трубы
		0.5	! позиция разделителя
		0.2	! времы распада разрыва
		2.0	! плотность слева
		0.0	! скорость слева
		2.0	! давление слева
		1.0	! плотность справа
		0.0	! скорость справа
		1.0	! давление справа
		1000	! количество ячеек
	Выходные данные: out.dat (координата, плотность, скорость, давления, внутренняя энергия)
	
	Теория изложена в статье: 
	E.F. Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics". Springer-Verlag, 1997          
*/
#include <stdio.h>
#include <math.h>

#define EPS			1.0e-6	// точность решения нелинейного уравнения
#define MAX_ITER	20		// количество итераций для решения нелинейного уравнения

/***** Константы из показателя адиабаты *****/
double	dgamma = 1.4, 
		g1 = (dgamma - 1.0)/(2.0*dgamma),
		g2 = (dgamma + 1.0)/(2.0*dgamma),
		g3 = 2.0*dgamma/(dgamma - 1.0),
		g4 = 2.0/(dgamma - 1.0),
		g5 = 2.0/(dgamma + 1.0),
		g6 = (dgamma - 1.0)/(dgamma + 1.0),
		g7 = (dgamma - 1.0)/2.0,
		g8 = dgamma - 1.0;

/***** Параметры ударной трубы *****/
	// Параметры слева
double	dl,  // плотность
		ul,  // скорость
		pl,  // давление
		cl;  // скорость звука

	// Параметры справа
double	dr,  // плотность
		ur,  // скорость
		pr,  // давление
		cr;  // скорость звука

/**************************************************
 Получение плотности, давления и скорости в точке
**************************************************/
void sample(double &pm, double &um, double &s, 
			double &d, double &u, double &p)
{
	double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

	if(s <= um)
	{
		// точка слева от контактного разрыва
		if(pm <= pl)
		{
			// левая волна разрежения
			shl = ul - cl;

			if(s <= shl)
			{
				// точка слева
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
					// точка слева от контактного разрыва
					d = dl * pow(pm/pl,1.0/dgamma);
					u = um;
					p = pm;
				}
				else
				{
					// точка в волне разрежения
					u = g5 * (cl + g7*ul + s);
					c = g5 * (cl + g7*(ul-s));
					d = dl*pow(c/cl,g4);
					p = pl*pow(c/cl,g3);
				}
			}
		}
		else
		{
			// левая ударная волна
			pml = pm/pl;
			sl  = ul - cl * sqrt(g2*pml + g1);

			if( s <= sl )
			{
				// точка слева от ударной волны
				d = dl;
				u = ul;
				p = pl;
			}
			else
			{
				// точка слева от контактного разрыва
				d = dl * (pml + g6)/(pml*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
	}
	else
	{
		// точка справа от контактного разрыва
		if( pm > pr )
		{
			// правая ударная волна
			pmr = pm/pr;
			sr  = ur + cr * sqrt(g2*pmr + g1);

			if( s >= sr )
			{
				// точка справа от ударной волны
				d = dr;
				u = ur;
				p = pr;
			}
			else
			{
				// точка справа от контактного разрыва
				d = dr * (pmr + g6)/(pmr*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
		else
		{
			// правая волна разрежения
			shr = ur + cr;

			if( s >= shr )
			{
				// точка справа от волны разрежения
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
					// точка справа от контактного разрыва
					d = dr * pow(pm/pr,1.0/dgamma);
					u = um;
					p = pm;
				}
				else
				{
					// точка в левой волне разрежения
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
		Решение для одной из частей
**************************************************/
void prefun(double &f, double &fd, double &p, 
			double &dk, double &pk, double &ck)
{
	double ak, bk, pratio, qrt; 
	
	if(p <= pk)
	{
		// волна разрежения
		pratio = p/pk;
		f	   = g4*ck*(pow(pratio,g1) - 1.0);
		fd     = (1.0/(dk*ck))*pow(pratio,-g2);

	}
	else
	{
		// ударная волна
		ak  = g5/dk;
		bk  = g6*pk;
		qrt = sqrt(ak/(bk + p));
		f   = (p - pk)*qrt;
		fd  = (1.0 - 0.5*(p - pk)/(bk + p)) * qrt;
	}

}

/**************************************************
		Вычисление начального приближения
**************************************************/
double guessp()
{
	double	cup, gel, ger, 
			pmax, pmin, ppv, pq, pm, 
			ptl, ptr, 
			qmax, quser, um;

	/*** Вычисление приближения давления из PVRS решения Римана ***/
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
			// две волны разрежения
            pq  = pow(pl/pr,g1);
            um  = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
            ptl = 1.0 + g7*(ul - um)/cl;
            ptr = 1.0 + g7*(um - ur)/cr;
            pm  = 0.5*(pl*pow(ptl,g3) + pr*pow(ptr,g3));
		}
		else
		{
			// две ударных волны
			gel = sqrt((g5/dl)/(g6*pl + ppv));
			ger = sqrt((g5/dr)/(g6*pr + ppv));
            pm  = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
		}

	}
	return pm;
}

/************************************************** 
Определение давления и скорости контактного разрыва			 
**************************************************/
void starpu(double &p, double &u)
{
	int i;
	
	double pstart,	// начальное приближение давления 
		   pold,	// предыдущее приближение давления
		   udiff,	// разность скоростей
		   change;	// разность давлений

	double	fl, fld, // функции давления
			fr, frd;	

	/*** Вычисление начального приближения ***/
	pstart = guessp();

	/*** Предыдущее приближение ***/
	pold = pstart;

	/*** Разность скоростей ***/
	udiff = ur - ul;

	/*** Метод Ньютона для определения давления ***/

	// разность между разными приближениями давлений
	change = 10.0 * EPS;

	for(i=0 ; i<MAX_ITER && change>EPS ; i++)
	{
		// решение для левой части
		prefun(fl,fld,pold,dl,pl,cl);

		// решение для правой части
		prefun(fr,frd,pold,dr,pr,cr);

		// очередное приближение давления
		p = pold - (fl + fr + udiff) / (fld + frd);

		// разность между разными приближениями давлений
		change = 2.0 * fabs((p-pold)/(p+pold));

		// если давление отрицательное, до обнуляем его
		if(p < 0.0) p =0.0;

		pold = p;
	}

	// определение скорости
	u = 0.5 * (ul + ur + fr - fl);
}

/************************************************** 
				Основная функция 
**************************************************/
int main()
{
	char string[128];	// строка в файле конфигурации ударной трубы

	double	diaph,	// позиция разделителя
			len,	// размер области [0;len]
			time,	// время распада разрыва
			pm,		// давление по сторонам контактного разрыва
			um,		// скорость контактного разрыва
			dx,		// шаг по пространству
			xpos,	// текущая точка области
			ds,		// плотность в точке области
			ps,		// давление в точке области
			us,		// скорость в точке области
			s;		// характеристика

	int numcells;	// количество ячеек
	int i;			// номер ячейки

	FILE *fout,		// файл результатов
		 *fin;		// файл с данными

	//*** задание параметров трубы ***
		fin = fopen("shock_tube.dat","r");
		fgets(string,sizeof(string),fin);
	// размер трубы [0;len]
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&len);
	// позиция разделителя
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&diaph);
	// время распада разрыва
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&time);
	// параметры слева
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&dl);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&ul);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&pl);

		cl	= sqrt(dgamma*pl/dl);
	// параметры справа
		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&dr);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&ur);

		fgets(string,sizeof(string),fin);
		sscanf(string,"%lf",&pr);

		cr	= sqrt(dgamma*pr/dr);
	// количество ячеек
		fgets(string,sizeof(string),fin);
		sscanf(string,"%d",&numcells);

	//********************************


	/*** Образование области вакуума ***/
	if( g4*(cl+cr) <= (ur-ul) )
	{
		printf("*** Generated vacuum ***\n");
		return -1;
	}

	/*** Нахождение точного решения давления и скорости контактного разрыва ***/
	starpu(pm,um); 

	/*** Определение шага (размера ячейки) ***/
	dx = len/double(numcells);

	/*** Вывод решения ***/
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
