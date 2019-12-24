# Sistema-Solar
Este programa de consola simula la dinámica del Sistema Solar (o en general de cualquier sistema de cuerpos sometidos únicamente a las fuerzas gravitatiorias que se ejercen entre ellos). Además, se pueden calcular la energía total, momento angular total y otras magnitudes de interés.

La simulación usa el algoritmo de Verlet (https://es.wikipedia.org/wiki/Integraci%C3%B3n_de_Verlet) y aunque el propósito era simular el Sistema Solar, se pueden añadir tantos cuerpos como se quieran y admite hacerlo en 1, 2, 3 o cualquier dimensión. Puede ser útil para simular cualquier sistema de cuerpos celestes o para dinámica molecular.

Los datos de cada planeta se cogen de un fichero (como por ejemplo se incluye el fichero Datos.dat) y los resultados pueden volcarse a un fichero o a un programa de visualización (que no se incluye en este proyecto). En este último caso, cabe reseñar que se han reescalado las ecuaciones ya que al tratar números muy grandes (masas en Kg) junto con números muy pequeños (la Constante Gravitatoria Universal) surgirían problemas con las aproximaciones. Concretamente se han reescalado de la siguiente forma:
    
    -Posiciones: r' = r/c ; con c = 1.496 × 10¹¹ m (distancia Tierra-Sol)
    -Tiempos: t' = t * [G*Ms / c³]^0.5 ; con Ms siendo la masa del Sol.
    -Masas: m' = m/Ms
    
Se debe tener en cuenta que los resultados que muestra la simulación están en estas unidades, así que es necesario hacer el proceso inverso para obtenerlos en unidades del Sistema Internacional.
