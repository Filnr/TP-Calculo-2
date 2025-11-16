# Sistema de Integração por Frações Parciais

Trabalho Prático desenvolvido para a disciplina de calculo 2
Professor(a): Cláudia Silva Tavares

Sistema computacional para resolução de integrais racionais utilizando o método de frações parciais.

## Índice
- [Sobre o Projeto](#sobre-o-projeto)
- [Requisitos](#requisitos)
- [Como Executar](#como-executar)
- [Modos de Uso](#modos-de-uso)
- [Tipos de Integrais Suportadas](#tipos-de-integrais-suportadas)
- [Exemplos de Uso](#exemplos-de-uso)
- [Estrutura do Código](#estrutura-do-código)

---

## Sobre o Projeto

Este sistema resolve integrais do tipo:

$$\int \frac{Ax + B}{\text{denominador}} \ dx$$

Onde o denominador pode ser:
- Produto de fatores lineares: `(x - x₁)(x - x₂)`
- Quadrático irredutível: `ax² + bx + c`
- Produto de quadráticos: `(x² + p₁)(x² + p₂)`

O programa decompõe automaticamente a fração em frações parciais e calcula a integral simbólica.

---

## Requisitos

- **Python 3.6+**
- Bibliotecas padrão: `math`, `re` (já incluídas no Python)
- **Nenhuma biblioteca externa necessária!**

---

## Como Executar

### 1. Clone o repositório
```bash
git clone https://github.com/Filnr/TP-Calculo-2.git
```

### 2. Execute o programa no terminal
```bash
python TP.py
```

### 3. Escolha o modo de operação
Ao executar, você verá o menu:
```
======================================================================
BEM-VINDO AO SISTEMA DE INTEGRAÇÃO POR FRAÇÕES PARCIAIS
======================================================================

Escolha uma opção:
1. Modo interativo (inserir valores manualmente)
2. Executar exemplos automáticos
3. Sair

Opção: 
```

---

## Modos de Uso

### **Modo 1: Interativo**
Permite inserir suas próprias integrais
Para expressar potência, utilize '^', como nos exemplos abaixo.

**Exemplo de sessão:**
```
--- NUMERADOR ---
Digite o numerador (ex: 3x + 5): x+2

--- DENOMINADOR ---
Digite o denominador (ex: x^2 - 3x + 2 ou (x-1)(x-2)): (x^2+1)(x^2+4)
```

### **Modo 2: Exemplos Automáticos**
Executa uma série de exemplos pré-configurados para demonstração.

---

## Tipos de Integrais Suportadas

### **Tipo 1: Fatores Lineares**
$$\int \frac{Ax + B}{(x - x_1)(x - x_2)} \, dx$$

**Entrada aceita:**
- Expandido: `x^2 - 3x + 2`
- Fatorado: `(x-1)(x-2)`
- Fatorado geral: `(3x-2)(x+5)`

**Decomposição:**
$$\frac{A_1}{x - x_1} + \frac{A_2}{x - x_2}$$

**Resultado:**
$$A_1 \ln|x - x_1| + A_2 \ln|x - x_2| + C$$

---

### **Tipo 2: Quadrático Irredutível**
$$\int \frac{Ax + B}{ax^2 + bx + c} \, dx$$

**Entrada aceita:**
- `x^2 + 2x + 5`
- `2x^2 + 3x + 7`

**Para raízes complexas (Δ < 0):**
$$\frac{C}{2a} \ln|ax^2 + bx + c| + D \arctan\left(\frac{x + h}{k}\right) + C$$

---

### **Tipo 3: Produto de Quadráticos**
$$\int \frac{Ax + B}{(x^2 + p_1)(x^2 + p_2)} \, dx$$

**Entrada aceita:**
- Expandido: `x^4 + 5x^2 + 4`
- **Fatorado: `(x^2+1)(x^2+4)`** 

**Decomposição:**
$$\frac{C_1x + D_1}{x^2 + p_1} + \frac{C_2x + D_2}{x^2 + p_2}$$

**Resultado:**
Combinação de logaritmos e arctangentes.

---

## Exemplos de Uso

### **Exemplo 1: Fatores Lineares Simples**
```
Numerador: 3x + 5
Denominador: (x-1)(x-2)

Resultado:
∫ f(x) dx = 8·ln|x - (1)| + -5·ln|x - (2)| + C
```

---

### **Exemplo 2: Fatores Lineares Gerais**
```
Numerador: 3x + 27
Denominador: (3x-2)(x+5)

Decomposição:
  Termo 1: 2.143 / (3x + -2)
  Termo 2: 0.8571 / (1x + 5)

Resultado:
∫ f(x) dx = 0.7143·ln|3x + -2| + 0.8571·ln|1x + 5| + C
```

---

### **Exemplo 3: Quadrático com Raízes Complexas**
```
Numerador: 2x + 3
Denominador: x^2 + 2x + 5

Resultado:
∫ f(x) dx = 1·ln|1x² + 2x + 5| + 0.5·arctan((x + 1) / 2) + C
```

---

### **Exemplo 4: Produto de Quadráticos** 
```
Numerador: x + 2
Denominador: (x^2+1)(x^2+4)

Decomposição:
  Termo 1: (0.3333x + 0.6667) / (x² + 1)
  Termo 2: (-0.3333x + 1.333) / (x² + 4)

Integração:
  ∫ Termo 1 dx = 0.1667·ln|1x² + 1| + 0.6667·arctan((x + 0) / 1)
  ∫ Termo 2 dx = -0.1667·ln|1x² + 4| + 0.6667·arctan((x + 0) / 2)

Resultado:
∫ f(x) dx = [combinação dos termos acima] + C
```

---

## Estrutura do Código

### **Classes Principais**

#### `Polinomio`
Representa polinômios através de coeficientes.
```python
p = Polinomio([2, -3, 1])  # Representa x² - 3x + 2
print(p.grau())  # Output: 2
```

#### `SistemaLinear`
Resolve sistemas lineares pelo método de Gauss.
```python
matriz = [[1, 1], [-2, -1]]
vetor = [3, -5]
solucao = SistemaLinear.resolver(matriz, vetor)
```

---

### **Funções Principais**

#### `parsePolinomio(expressao)`
Converte string em coeficientes ou detecta forma fatorada.
```python
parsePolinomio("3x + 5")        # → [5, 3]
parsePolinomio("(x-1)(x-2)")   # → {'fatorado': True, ...}
```

#### `verificaFracao(numerador, denominador)`
Valida se a fração é apropriada para o método.
- Verifica fração própria
- Checa graus válidos
- Retorna mensagens de erro claras

#### `identificaTipoFatoracao(denominador)`
Identifica automaticamente o tipo de fatoração:
- `linear`: duas raízes reais
- `quadratico_complexo`: raízes complexas
- `misto`: produto de quadráticos
- `linear_fatorado`: mantém forma (ax+b)(cx+d)
- `misto_fatorado`: mantém forma (x²+p₁)(x²+p₂)

#### `decompoeEmFracoesParciais(...)`
Decompõe a fração usando sistemas lineares.

#### `integraCadaTermo(termo)`
Integra cada termo da decomposição.

#### `calculaIntegral(numerador, denominador)`
**Função coordenadora principal** - executa todo o pipeline.

---

## Formato de Entrada

### **Notações Aceitas**

| Tipo | Exemplos Válidos |
|------|------------------|
| Linear | `3x + 5`, `-2x + 7`, `x - 3` |
| Quadrático | `x^2 + 2x + 1`, `2x^2 - 3x + 5` |
| Quártico | `x^4 + 5x^2 + 4` |
| Fatorado Linear | `(x-1)(x-2)`, `(3x-2)(x+5)` |
| Fatorado Quadrático | `(x^2+1)(x^2+4)` |

### **Operadores**
- Potência: `x^2` ou `x**2`
- Multiplicação implícita: `3x` (não precisa `3*x`)
- Espaços opcionais

---

## Limitações

### **Não Suportado:**
- Frações impróprias (grau numerador ≥ grau denominador)
- Numeradores com grau > 1
- Denominadores com grau ≠ 2 ou 4
- Fatores repetidos (ex: `(x-1)²`)

### **Solução para Frações Impróprias:**
Use divisão polinomial primeiro, depois integre o quociente e o resto separadamente.

---

## Troubleshooting

### **Erro: "Fração imprópria"**
**Causa:** Grau do numerador ≥ grau do denominador.
**Solução:** Divida os polinômios antes.

### **Erro: "Grau do denominador não suportado"**
**Causa:** Denominador com grau diferente de 2 ou 4.
**Solução:** Fatore ou simplifique o denominador.

### **Decomposição vazia**
**Causa:** Tipo de fatoração não reconhecido.
**Solução:** Verifique o formato de entrada, use forma fatorada explícita.

---

## Teoria Matemática

### **Método das Frações Parciais**

Para integrar $\frac{P(x)}{Q(x)}$ onde $\text{grau}(P) < \text{grau}(Q)$:

1. **Fatore** $Q(x)$ em fatores irredutíveis
2. **Decomponha** em soma de frações simples
3. **Resolva** sistema linear para encontrar coeficientes
4. **Integre** cada fração simples

### **Sistemas de Equações**

O código resolve automaticamente sistemas lineares usando:
- **Eliminação de Gauss** com pivoteamento parcial
- **Substituição reversa**

---

## Autores

Filipe Nery Rabelo
Lucca Sander FrissoLucca Sander Frisso 
---

## Licença

Projeto acadêmico - Uso educacional livre.

---

## Suporte

Para dúvidas ou problemas:
1. Verifique os exemplos na seção [Exemplos de Uso](#exemplos-de-uso)
2. Execute o modo 2 (exemplos automáticos) para ver casos funcionando
3. Consulte a seção [Troubleshooting](#troubleshooting)
