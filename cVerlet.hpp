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


/*Esta clase es la base de la simulación, es la que hace todos los cálculos y contiene el
algoritmo de Verlet propiamente dicho.*/

#include <iostream>
#include <string>

class cVerlet
{
	public:
		double h; //h paso de la simulación.
		int centrarEnObj; //En qué objeto se centra el sistema de coordenadas (si vale negativo no se centra en ninguno).

		
		cVerlet(int nobj, int ndim); 
		~cVerlet();

		void cargar(std::string fichero); //Carga los datos de un fichero.
		void reescalar(double dm, double dr, double dt, double nuevaG); //Multiplica cada magnitud por las constantes dadas.
		
		void mostrarPosiciones (); //Muestra las posiciones actuales.
		double mostrarEnergia (); //Muestra la energía total actual.
		double mostrarMomenAng (); //Muestra el momento angular total actual.
		void mostrarPeriodo (int obj, double Ax); //Muestra el periodo del objeto dado.

		void arrancar(int iter); //Arranca la simulación y calcula el número de iteraciones dada. Muestra las posiciones automáticamente
		void calcularPaso (); //Calcula un único paso pero no muestra nada.


	private:
		int n, dim;
		double G;
		//Definimos las variables básicas. El primer corchete es el cuerpo, y el segundo es la componente.
		double **r, **v, **a, *m; //Son matrices dinámicas.
		void calcPos (); //Calcula las nuevas posiciones de cada cuerpo.
		void calcVel (); //Calcula las nuevas velocidades de cada cuerpo.
		void calcAcelG ();//Calcula las aceleraciones por gravedad que experimentan todos los cuerpos.
};
