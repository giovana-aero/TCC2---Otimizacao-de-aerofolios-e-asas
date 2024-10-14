# TCC2---Otimizacao-de-aerofolios-e-asas

Compilado dos scripts utilizados no meu trabalho de conclusão de curso do curso de Engenharia Aeroespacial na Universidade de Brasília. A proposta foi o desenvolvimento de algoritmos de otimização de aerofólios e asas utilizando aloritmos genéticos escritos em Octave e Python.

Os scripts .m foram escritos no Octave. Não foram feitos testes no MATLAB, mas devem funcionar sem problema.

Os executáveis relevantes do XFOIL, do APAME e do VSPAERO já estão inclusos nestas pastas, mas, para fim de documentação, eles estão disponíveis em, respectivamente, https://web.mit.edu/drela/Public/web/xfoil/, http://www.3dpanelmethod.com/download.html e http://openvsp.org/download.php

Há alguns aspectos destes códigos que eu não pude testar completamente, portanto, erros podem ocorrer.

<< Autora: Giovana Fernandes >>


A ser alterado:
* Adicionar uma forma de mudar o número de Mach nos algoritmos 2D (XFOIL)
* Adicionar uma forma de alterar o ponto de referência para o cálculo do CM - no XFOIL isso é feito com o comando XYCM; no APAME é preciso transladar a geometria em relação à origem; no VSPAERO há comandos específicos na API para isso
