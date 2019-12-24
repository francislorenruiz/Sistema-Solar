/*Copyright 2017 Francisco Lorente Ruiz.

    This file is part of Sistema Solar.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>. 
*/


#include "cVerlet.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <cmath>

cVerlet::cVerlet(int nobj, int ndim)
{
	int i, j;

	n=nobj;
	dim=ndim;
	centrarEnObj=-1;
	h=0.01;
	G=6.67E-11; //Constante gravitatoria universal en unidades SI.


	//Iniciamos las matrices dinámicas.
	r = new double *[n];
	v = new double *[n];
	a = new double *[n];
	m = new double [n];
	for(i=0; i<n; i++)
	{
		r[i] = new double[dim];
		v[i] = new double[dim];
		a[i] = new double[dim];
	}
}

cVerlet::~cVerlet()
{
	int i;

	//Destruimos las matrices dinámicas
	for(i=0; i<n; i++)
	{
		delete[] r[i];
		delete[] v[i];
		delete[] a[i];
	}

	delete[] r;
	delete[] v;
	delete[] a;
	delete[] m;
}

void cVerlet::cargar(std::string fichero)
{
	int i, j;
	std::ifstream fich;

	//Leemos los datos de un fichero
	fich.open(fichero.c_str());
	fich.ignore(std::numeric_limits<int>::max(), '\n'); //Esta orden ignora la primera línea del fichero de datos. Lo reservo para poner que es cada número y que esté más claro.

	//Ahora sí, lee los datos.
	for(i=0; i<n; i++)
	{
		fich>>m[i];
		for(j=0; j<dim; j++) fich>>r[i][j];
		for(j=0; j<dim; j++) fich>>v[i][j];
	}

	fich.close();

	//Si se quiere centrar el sistema de referencia en un objeto particular se desplazan todos los objetos para que el deseado esté en el (0,0,0).
	if(centrarEnObj>=0)
	{	
		for(i=0; i<n; i++)
			if(i!=centrarEnObj)
				for(j=0; j<dim; j++)
					r[i][j]-=r[centrarEnObj][j];
		for(j=0; j<dim; j++) r[centrarEnObj][j]=0;
	}

	return;
}

void cVerlet::reescalar(double dm, double dr, double dt, double nuevaG)
{
	int i, j;
	double aux=dr/dt;

	for(i=0; i<n; i++)
	{
		m[i]*=dm;
		for(j=0; j<dim; j++)
		{
			r[i][j]*=dr;
			v[i][j]*=aux;
		}
	}

	G=nuevaG;

	return;
}

void cVerlet::mostrarPosiciones ()
{
	int i, j;	

	for(i=0; i<n; i++)
	{
		for(j=0; j<dim; j++)
			std::cout<<r[i][j]<<" ";
		std::cout<<std::endl;
	}
	std::cout<<std::endl;

	return;
}

double cVerlet::mostrarEnergia ()
{
	int i, j, k;
	double dist, energia=0;

	for(i=0; i<n; i++) //Calculamos la energía cinética
		for(j=0; j<dim; j++)
			energia+=1.0/2*m[i]*v[i][j]*v[i][j];

	for(i=0; i<n-1; i++) //Calculamos la energía potencial
	{
		for(j=i+1; j<n; j++)
		{
			dist=0; //Calculamos la distancia entre los planetas
			for(k=0; k<dim; k++) dist+=(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
			dist=sqrt(dist);

			energia-=G*m[i]*m[j]/dist;
		}
	}

	std::cout <<energia <<std::endl;

	return energia;
}

double cVerlet::mostrarMomenAng ()
{
	int i, j;
	double *p, L=0;

	p = new double[dim];

	for(i=0; i<n; i++)
	{
		for(j=0; j<dim; j++) p[j]=m[i]*v[i][j];
		L+=r[i][0]*p[1]-r[i][1]*p[0];
	}

	std::cout <<L <<std::endl;

	delete[] p;
	
	return L;
}

void cVerlet::mostrarPeriodo (int obj, double Ax)
{
	int i;
	double periodo=0.0, *r0;
	bool acabar=false, salirDelIntervalo=false;
	
	r0 = new double[dim];

	for(i=0; i<dim; i++) r0[i]= r[obj][i];

	do{
		calcularPaso();
		for(i=0; i<dim; i++)
		{
			acabar=false;
			if(r0[i]+Ax>=r[obj][i] && r0[i]-Ax<=r[obj][i]) 
			{			
				if(salirDelIntervalo==true)acabar=true;
			}
			else
			{	
				salirDelIntervalo=true;
				acabar=false;
				break;
			}
		}
		if(acabar) periodo*=h;
		else periodo++;
	}while(!acabar);

	std::cout<<periodo<<std::endl;

	delete[] r0;

	return;
}

void cVerlet::arrancar(int iter)
{
	int i, j;	
	
	mostrarPosiciones();
	calcAcelG();
	for(i=0; i<iter; i++)
	{
		calcPos ();
		calcVel ();
		calcAcelG();
		calcVel ();

		mostrarPosiciones();
	}

	return;
}

void cVerlet::calcularPaso()
{
	calcAcelG();

	calcPos ();
	calcVel ();
	calcAcelG();
	calcVel ();

	return;
}


void cVerlet::calcPos ()
{
	int i, j;

	for(i=0; i<n; i++)
		for(j=0; j<dim; j++) 
			r[i][j]+= h*v[i][j] + h*h*a[i][j]/2.0;
	
	if(centrarEnObj>=0)
	{	
		for(i=0; i<n; i++)
			if(i!=centrarEnObj)
				for(j=0; j<dim; j++)
					r[i][j]-=r[centrarEnObj][j];
		for(j=0; j<dim; j++) r[centrarEnObj][j]=0;
	}

	return;
}

void cVerlet::calcVel ()
{
	int i, j;

	for(i=0; i<n; i++)
		for(j=0; j<dim; j++) 
			v[i][j]+= h*a[i][j]/2.0;

	return;
}

void cVerlet::calcAcelG ()
{
	int i, j, k;
	double dist, aux;

	//Hacemos cero todas las aceleraciones.
	for(i=0; i<n; i++) 
		for(j=0; j<dim; j++)
			a[i][j]=0;	
	
	for(i=0; i<n-1; i++) //Para cada planeta se calcula su aceleración.
	{
		for(j=i+1; j<n; j++) //Nos movemos entre los planetas
		{
			dist=0; //Calculamos la distancia entre los planetas y elevamos al cubo
			for(k=0; k<dim; k++) dist+=(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
			dist=pow(dist, 1.5);
			
			//Nos movemos entre las componentes	
			for(k=0; k<dim; k++) 
			{
				aux=G*(r[i][k]-r[j][k]) / dist;
				a[i][k]-=m[j]*aux;
				a[j][k]+=m[i]*aux; //La fuerza que ejerce una partícula sobre otra es igual pero de signo contrario que la que ejerce la otra sobre la una.
			}
		}
	}

	return;
}

