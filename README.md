# HW3_P3

Mallas quad9 utilizadas:

![image](https://user-images.githubusercontent.com/81662690/131737423-4ad9d855-1d3f-487d-aa91-d5ad17cbee79.png)
![image](https://user-images.githubusercontent.com/81662690/131737444-d382a1ba-eb69-41ad-9233-5657e3d7d375.png)
![image](https://user-images.githubusercontent.com/81662690/131737470-920f312d-ab20-402e-a17e-6cc21dbd8f93.png)

Para las deformaciones se amplificó por 5e7

PARTE A:

![image](https://user-images.githubusercontent.com/81662690/131735837-e43e5262-513c-4879-b0fd-f1c83e518a4e.png)

Para esta parte de la tarea se utilizaron los datos obtenidos para le quad4 de la parte anterior y los obtenidos al desarrollar esta parte, se pueden comprobar las tensiones máximas más adelante. El valor "H" es la distancia mínima entre nodos, esta varía de un modelo a otro siendo que se utilizó casi exáctamente la misma cantidad de nodos pero cambia por la distribución de estos. A continuación se muestra el gráfico:

![image](https://user-images.githubusercontent.com/81662690/131736823-9c092c75-6ef1-4c82-a73b-eb0e822693ae.png)

Se puede ver claramente que el modelo quad9 converge antes por lo que se puede usar una malla menos detallada utilizando un menor recurso computacional. La malla media del quad9 alcanza el máximo (el cual se mantiene en la más detallada) mientras que en la malla media del quad4 aun está lejos de converger siendo que tienen un parámetro H similar y una cantidad de nodos casi idéntica.

Para efectos de asegurar convergencia y optimización del recurso computacional preferiría utilizar el modelo quad9 ya que su implementación es similar y para este caso entrega mejores resultados.


PARTE B:

![image](https://user-images.githubusercontent.com/81662690/131737332-7be3dced-c0b9-447b-9763-343372b2b4c1.png)

Esta parte se realizó solo para el modelo quad9 ya que eso acordamos con el profesor, aun que se implementa de la misma forma para el modelo quad4 (cambiando la matriz de elasticidad). Se muentran las diferencias de desplazamientos para la malla b (refimamiento medio) ya que en la parte anterior se comprobó que esta malla ya alcanzaba la tensión máxima.

Malla c:
Deformaciones
alpha = 1
![image](https://user-images.githubusercontent.com/81662690/131738939-9274263b-088c-4343-95a7-5aabc81cfc6e.png)

alpha = 2
![image](https://user-images.githubusercontent.com/81662690/131739171-a60530c1-a277-408c-aa32-18357cc68217.png)

alpha = 4
![image](https://user-images.githubusercontent.com/81662690/131739214-dcc07843-6716-42bd-bad3-a710583bc62b.png)

Tensiones
alpha = 1
![image](https://user-images.githubusercontent.com/81662690/131739482-c63c8a92-5394-4ba6-9315-06a6f9bde05f.png)

alpha = 2
![image](https://user-images.githubusercontent.com/81662690/131739575-34588d80-be03-418c-9bc7-6107f31dec8d.png)

alpha = 4
![image](https://user-images.githubusercontent.com/81662690/131739623-206d4d64-29e6-406b-995b-f57007f62eba.png)

Malla b:
Deformaciones
alpha = 1
![image](https://user-images.githubusercontent.com/81662690/131740210-2bdc02f9-da30-48de-a8eb-b9503cd3fd10.png)

alpha = 2
![image](https://user-images.githubusercontent.com/81662690/131740257-27a51c5e-549b-410a-b5fc-50efa4da002b.png)

alpha = 4
![image](https://user-images.githubusercontent.com/81662690/131740314-07c848cd-94ae-4e64-aab8-738a033b01a3.png)

Tensiones
alpha = 1
![image](https://user-images.githubusercontent.com/81662690/131740518-3fa51462-b17c-4221-ab7c-56c0b478dbd7.png)

alpha = 2
![image](https://user-images.githubusercontent.com/81662690/131740630-caa832ed-c434-4d64-8459-fce314525954.png)

alpha = 4
![image](https://user-images.githubusercontent.com/81662690/131740692-f7e6ec98-2da4-4000-9bfc-edc4e6e5a765.png)

En los 3 casos se puede apreciar que aumenta levemente la deformación en la dirección longitudinal pero que la deformación cercana al orificio en la dirección "y" es mayor al aumentar el factor alpha. Además, al aumentar el factor alpha se puede apreciar que las tensiones máximas aumentan pero que además disminuye el área de influencia por lo que se concentran en un área menor, esto podría ser algo muy útil para aumentar en espesor en un área mas reducida y así optimizar el uso de material.

No coloqué la malla más detallada porque la conclusión sería la misma y es poca la diferencia con respecto a la malla intermedia.
