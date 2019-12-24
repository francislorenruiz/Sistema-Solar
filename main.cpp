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


#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>
#include "cVerlet.hpp"

using namespace std;

int main()
{
	int i;
	double G, c, ms;
	cVerlet simulacion(9, 2); //Creamos la clase y le damos el número de planetas y de dimensiones con las que trabajará.
	ofstream fichero;

	G=6.67E-11;
	c=1.496E11;
	ms=1.99E30;

	simulacion.h = 0.01;
	simulacion.cargar("Datos.dat");
	simulacion.reescalar(1.0/ms, 1.0/c, sqrt(G*ms/pow(c, 3)), 1.0);
	simulacion.centrarEnObj = -1; //Se centra en ningún objeto.

	simulacion.mostrarEnergia();
	for(i=0; i<100000; i++) 
	{
		simulacion.calcularPaso();
		simulacion.mostrarPosiciones();
	}

	return 0;
}
