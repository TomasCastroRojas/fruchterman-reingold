# Fruchterman-Reingold
Implementación del algoritmo de Fruchterman-Reingold para la visualización de grafos

## Uso
Para ejecutar escribir en una terminal
```code
python main.py grafo.txt
```
donde el archivo `texto.txt` debe tener el siguiente formato
```code
n # numero que indica cantidad de verticies
id_1
id_2
...
id_n
id_1 id_3 # Esto representa una arista desde el vertice id_1 al vertice id_3
id_n id_4
...
```
Las aristas no deben repetirse. Se pueden encontrar varios ejemplos en el directorio `grafos`.
Los distintos parametros del parametro influyen en como varia la visualizacion paso a paso. No hay una configuración ideal.
Los valores por defecto sirven para mostrar una animación fluida. Se alienta probar variaciones y ver cómo afecta al proceso.
### Dependencias
Este proyecto requiere `matplotlib`, `argparse` y `numpy`.