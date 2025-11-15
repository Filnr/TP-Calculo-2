"""
Trabalho Prático - Integração por Frações Parciais
Sistema interativo para resolução de integrais racionais
Implementação sem SymPy e NumPy
"""

import math
import re


# ==============================================================================
# CLASSE AUXILIAR: POLINÔMIO
# ==============================================================================

class Polinomio:
    """Representa um polinômio através de seus coeficientes"""
    
    def __init__(self, coefs):
        """coefs[i] é o coeficiente de x^i"""
        self.coefs = list(coefs)
        self._limpar()
    
    def _limpar(self):
        """Remove zeros à direita"""
        while len(self.coefs) > 1 and abs(self.coefs[-1]) < 1e-10:
            self.coefs.pop()
    
    def grau(self):
        """Retorna o grau do polinômio"""
        return len(self.coefs) - 1
    
    def __str__(self):
        if not self.coefs:
            return "0"
        
        termos = []
        for i in range(len(self.coefs) - 1, -1, -1):
            c = self.coefs[i]
            if abs(c) < 1e-10:
                continue
            
            sinal = "+" if c > 0 and termos else ""
            if c < 0:
                sinal = "-"
            
            valor = abs(c)
            
            if i == 0:
                termos.append(f"{sinal}{valor:.4g}")
            elif i == 1:
                if abs(valor - 1) < 1e-10:
                    termos.append(f"{sinal}x")
                else:
                    termos.append(f"{sinal}{valor:.4g}x")
            else:
                if abs(valor - 1) < 1e-10:
                    termos.append(f"{sinal}x^{i}")
                else:
                    termos.append(f"{sinal}{valor:.4g}x^{i}")
        
        return "".join(termos) if termos else "0"


# ==============================================================================
# CLASSE AUXILIAR: SISTEMA LINEAR
# ==============================================================================

class SistemaLinear:
    """Resolve sistemas lineares pelo método de Gauss"""
    
    @staticmethod
    def resolver(matriz, vetor):
        """
        Resolve Ax = b pelo método de eliminação de Gauss
        matriz: matriz A (lista de listas)
        vetor: vetor b (lista)
        Retorna: vetor solução x
        """
        n = len(vetor)
        # Cria matriz aumentada [A|b]
        A = [linha[:] + [vetor[i]] for i, linha in enumerate(matriz)]
        
        # Eliminação de Gauss com pivoteamento parcial
        for i in range(n):
            # Pivoteamento
            max_linha = i
            for k in range(i + 1, n):
                if abs(A[k][i]) > abs(A[max_linha][i]):
                    max_linha = k
            A[i], A[max_linha] = A[max_linha], A[i]
            
            # Eliminação
            for k in range(i + 1, n):
                if abs(A[i][i]) < 1e-10:
                    continue
                fator = A[k][i] / A[i][i]
                for j in range(i, n + 1):
                    A[k][j] -= fator * A[i][j]
        
        # Substituição reversa
        x = [0] * n
        for i in range(n - 1, -1, -1):
            if abs(A[i][i]) < 1e-10:
                x[i] = 0
            else:
                x[i] = A[i][n]
                for j in range(i + 1, n):
                    x[i] -= A[i][j] * x[j]
                x[i] /= A[i][i]
        
        return x


# ==============================================================================
# FUNÇÃO DE PARSING
# ==============================================================================

def parsePolinomio(expressao):
    """
    Converte uma string de polinômio em lista de coeficientes.
    
    Exemplos:
        "3x + 5" → [5, 3] (coefs de x^0, x^1)
        "x^2 - 3x + 2" → [2, -3, 1]
        "10x^2 - 20" → [-20, 0, 10]
    
    Retorna:
        lista de coeficientes [c0, c1, c2, ...] onde ci é coef de x^i
        OU
        {'fatorado': True, 'fatores': [...]} se estiver em forma fatorada
    """
    # Remove espaços
    expressao = expressao.replace(" ", "")
    
    # Substitui formas comuns
    expressao = expressao.replace("**", "^")
    
    # Verifica se há produto de fatores lineares (forma fatorada)
    if expressao.count('(') >= 2 and expressao.count(')') >= 2:
        return extrairFatores(expressao)
    
    # Adiciona sinal + no início se não houver
    if not expressao:
        return [0]
    if expressao[0] not in ['+', '-']:
        expressao = '+' + expressao
    
    # Encontra todos os termos usando regex melhorado
    # Padrão: sinal + coeficiente (opcional) + x (opcional) + ^expoente (opcional)
    padrao = r'([+-])\s*(\d*\.?\d*)\s*\*?\s*([xX])?\s*\^?\s*(\d*)'
    termos = re.findall(padrao, expressao)
    
    # Dicionário para armazenar coeficientes por grau
    coefs_dict = {}
    
    for sinal, coef, var, exp in termos:
        if not coef and not var:
            continue
        
        # Determina o sinal
        multiplicador = 1 if sinal == '+' else -1
        
        # Determina o coeficiente
        if var and not coef:
            coef_valor = 1.0
        elif not var and coef:
            coef_valor = float(coef)
        elif coef:
            coef_valor = float(coef)
        else:
            continue
        
        coef_valor *= multiplicador
        
        # Determina o expoente
        if not var:
            expoente = 0
        elif not exp:
            expoente = 1
        else:
            expoente = int(exp)
        
        # Acumula coeficientes do mesmo grau
        if expoente in coefs_dict:
            coefs_dict[expoente] += coef_valor
        else:
            coefs_dict[expoente] = coef_valor
    
    # Converte para lista
    if not coefs_dict:
        return [0]
    
    grau_max = max(coefs_dict.keys())
    coeficientes = [0] * (grau_max + 1)
    
    for grau, coef in coefs_dict.items():
        coeficientes[grau] = coef
    
    return coeficientes


def extrairFatores(expressao):
    """
    Extrai fatores individuais sem expandir.
    Exemplos:
        "(3x-2)(x+5)" → {'fatorado': True, 'fatores': [[coefs1], [coefs2]]}
    """
    # Extrai os fatores entre parênteses
    fatores_str = re.findall(r'\([^)]+\)', expressao)
    
    if not fatores_str:
        return parsePolinomio(expressao)
    
    # Parseia cada fator
    fatores = []
    for fator_str in fatores_str:
        # Remove parênteses
        fator_limpo = fator_str[1:-1]
        coefs = parsePolinomioSimples(fator_limpo)
        fatores.append(coefs)
    
    return {'fatorado': True, 'fatores': fatores}


def parsePolinomioSimples(expressao):
    """
    Versão simplificada do parsePolinomio para fatores individuais.
    """
    expressao = expressao.replace(" ", "").replace("**", "^")
    
    if not expressao:
        return [0]
    if expressao[0] not in ['+', '-']:
        expressao = '+' + expressao
    
    padrao = r'([+-])\s*(\d*\.?\d*)\s*\*?\s*([xX])?\s*\^?\s*(\d*)'
    termos = re.findall(padrao, expressao)
    
    coefs_dict = {}
    
    for sinal, coef, var, exp in termos:
        if not coef and not var:
            continue
        
        multiplicador = 1 if sinal == '+' else -1
        
        if var and not coef:
            coef_valor = 1.0
        elif not var and coef:
            coef_valor = float(coef)
        elif coef:
            coef_valor = float(coef)
        else:
            continue
        
        coef_valor *= multiplicador
        
        if not var:
            expoente = 0
        elif not exp:
            expoente = 1
        else:
            expoente = int(exp)
        
        if expoente in coefs_dict:
            coefs_dict[expoente] += coef_valor
        else:
            coefs_dict[expoente] = coef_valor
    
    if not coefs_dict:
        return [0]
    
    grau_max = max(coefs_dict.keys())
    coeficientes = [0] * (grau_max + 1)
    
    for grau, coef in coefs_dict.items():
        coeficientes[grau] = coef
    
    return coeficientes


def multiplicarPolinomios(p1, p2):
    """
    Multiplica dois polinômios representados por suas listas de coeficientes.
    """
    grau_resultado = len(p1) + len(p2) - 2
    resultado = [0] * (grau_resultado + 1)
    
    for i in range(len(p1)):
        for j in range(len(p2)):
            resultado[i + j] += p1[i] * p2[j]
    
    return resultado


# ==============================================================================
# FUNÇÕES PRINCIPAIS
# ==============================================================================

def verificaFracao(numerador, denominador):
    """
    Verifica se a fração é válida para integração por frações parciais.
    
    Parâmetros:
        numerador: lista de coeficientes OU dicionário com forma fatorada
        denominador: lista de coeficientes OU dicionário com forma fatorada
    
    Retorna:
        (valido, mensagem, denominador_expandido)
    """
    # Se denominador está fatorado, expande para verificação
    den_expandido = denominador
    if isinstance(denominador, dict) and denominador.get('fatorado'):
        fatores = denominador['fatores']
        den_expandido = fatores[0]
        for i in range(1, len(fatores)):
            den_expandido = multiplicarPolinomios(den_expandido, fatores[i])
    
    # Se numerador está fatorado, expande
    num_expandido = numerador
    if isinstance(numerador, dict) and numerador.get('fatorado'):
        fatores = numerador['fatores']
        num_expandido = fatores[0]
        for i in range(1, len(fatores)):
            num_expandido = multiplicarPolinomios(num_expandido, fatores[i])
    
    num = Polinomio(num_expandido)
    den = Polinomio(den_expandido)
    
    # Verifica se denominador não é zero
    if den.grau() == 0 and abs(den.coefs[0]) < 1e-10:
        return False, "Erro: Denominador não pode ser zero", den_expandido
    
    # Verifica se é fração própria (grau numerador < grau denominador)
    if num.grau() >= den.grau():
        return False, f"Erro: Fração imprópria (grau numerador {num.grau()} >= grau denominador {den.grau()}). Use divisão polinomial primeiro.", den_expandido
    
    # Verifica grau do denominador (deve ser 2 ou 4 para os casos propostos)
    if den.grau() not in [2, 4]:
        return False, f"Erro: Grau do denominador ({den.grau()}) não suportado. Deve ser 2 ou 4.", den_expandido
    
    # Verifica grau do numerador (deve ser no máximo 1)
    if num.grau() > 1:
        return False, f"Erro: Grau do numerador ({num.grau()}) deve ser no máximo 1 (forma Ax + B).", den_expandido
    
    return True, f"Fração válida: ({num}) / ({den})", den_expandido


def identificaTipoFatoracao(denominador):
    """
    Identifica o tipo de fatoração do denominador.
    
    Parâmetros:
        denominador: lista de coeficientes OU dicionário com forma fatorada
    
    Retorna:
        (tipo, info)
    """
    # Se já está fatorado, identifica pelos fatores
    if isinstance(denominador, dict) and denominador.get('fatorado'):
        fatores = denominador['fatores']
        
        # Verifica se são todos lineares (grau 1)
        todos_lineares = all(len(f) == 2 for f in fatores)
        
        if todos_lineares and len(fatores) == 2:
            # (ax + b)(cx + d) - fatores lineares
            # Extrai as raízes de cada fator ax + b = 0 → x = -b/a
            raizes = []
            fatores_str = []
            for fator in fatores:
                b, a = fator[0], fator[1]
                if abs(a) > 1e-10:
                    raiz = -b / a
                    raizes.append(raiz)
                    fatores_str.append(f"({a:.4g}x + {b:.4g})")
            
            return 'linear_fatorado', {
                'descricao': 'Produto de fatores lineares (mantido fatorado)',
                'raizes': raizes,
                'fatores': fatores,
                'fatores_str': " × ".join(fatores_str)
            }
        
        # Se tem fatores quadráticos
        todos_quadraticos = all(len(f) == 3 for f in fatores)
        if todos_quadraticos and len(fatores) == 2:
            return 'misto_fatorado', {
                'descricao': 'Produto de fatores quadráticos (mantido fatorado)',
                'fatores': fatores
            }
    
    # Caso contrário, usa a lógica expandida original
    den = Polinomio(denominador) if isinstance(denominador, list) else Polinomio(denominador)
    grau = den.grau()
    
    if grau == 2:
        # ax² + bx + c
        c, b, a = den.coefs[0], den.coefs[1] if len(den.coefs) > 1 else 0, den.coefs[2]
        delta = b**2 - 4*a*c
        
        if delta > 0:
            # Duas raízes reais distintas
            x1 = (-b + math.sqrt(delta)) / (2*a)
            x2 = (-b - math.sqrt(delta)) / (2*a)
            return 'linear', {
                'descricao': 'Produto de fatores lineares (x - x₁)(x - x₂)',
                'raizes': [x1, x2],
                'fatores': f"({a:.4g})(x - {x1:.4g})(x - {x2:.4g})"
            }
        elif abs(delta) < 1e-10:
            # Raiz dupla
            x1 = -b / (2*a)
            return 'linear_dupla', {
                'descricao': 'Fator linear repetido (x - x₁)²',
                'raizes': [x1],
                'fatores': f"({a:.4g})(x - {x1:.4g})²"
            }
        else:
            # Raízes complexas
            return 'quadratico_complexo', {
                'descricao': 'Quadrático irredutível (raízes complexas)',
                'coeficientes': (a, b, c),
                'fatores': f"{a:.4g}x² + {b:.4g}x + {c:.4g}"
            }
    
    elif grau == 4:
        # Verifica se é produto de dois quadráticos (x² + x₁)(x² + x₂)
        coefs = den.coefs
        
        # Verifica se coeficientes de x³ e x são zero
        if len(coefs) > 1 and abs(coefs[1]) < 1e-10 and len(coefs) > 3 and abs(coefs[3]) < 1e-10:
            # É da forma ax⁴ + bx² + c
            c, b, a = coefs[0], coefs[2], coefs[4]
            
            # Resolve equação do segundo grau em x²: at² + bt + c = 0 onde t = x²
            delta_t = b**2 - 4*a*c
            
            if delta_t >= 0:
                t1 = (-b + math.sqrt(delta_t)) / (2*a)
                t2 = (-b - math.sqrt(delta_t)) / (2*a)
                
                if t1 > 0 and t2 > 0:
                    return 'misto', {
                        'descricao': 'Produto de quadráticos (x² + x₁)(x² + x₂)',
                        'valores': [t1, t2],
                        'fatores': f"(x² + {t1:.4g})(x² + {t2:.4g})"
                    }
        
        return 'complexo', {
            'descricao': 'Fatoração complexa de grau 4',
            'fatores': str(den)
        }
    
    return 'desconhecido', {'descricao': 'Tipo não identificado'}


def decompoeEmFracoesParciais(numerador, denominador, tipo_fatoracao, info_fatoracao):
    """
    Decompõe a fração em frações parciais.
    
    Retorna:
        lista de tuplas (coeficientes, denominador_info)
    """
    # Expande numerador se necessário
    num_expandido = numerador
    if isinstance(numerador, dict) and numerador.get('fatorado'):
        fatores = numerador['fatores']
        num_expandido = fatores[0]
        for i in range(1, len(fatores)):
            num_expandido = multiplicarPolinomios(num_expandido, fatores[i])
    
    num = Polinomio(num_expandido)
    A = num.coefs[1] if len(num.coefs) > 1 else 0
    B = num.coefs[0]
    
    decomposicao = []
    
    print(f"DEBUG: A={A}, B={B}, tipo={tipo_fatoracao}")
    
    if tipo_fatoracao == 'linear_fatorado':
        # Trabalha diretamente com os fatores (ax + b)(cx + d)
        # (Ax + B) / [(a₁x + b₁)(a₂x + b₂)] = A₁/(a₁x + b₁) + A₂/(a₂x + b₂)
        
        fatores = info_fatoracao['fatores']
        raizes = info_fatoracao['raizes']
        
        # Para (Ax + B) / [(a₁x + b₁)(a₂x + b₂)]
        # Usando método dos coeficientes indeterminados:
        # Ax + B = A₁(a₂x + b₂) + A₂(a₁x + b₁)
        
        # Coeficientes dos fatores
        b1, a1 = fatores[0][0], fatores[0][1]
        b2, a2 = fatores[1][0], fatores[1][1]
        
        # Sistema:
        # A₁·a₂ + A₂·a₁ = A
        # A₁·b₂ + A₂·b₁ = B
        
        matriz = [[a2, a1], [b2, b1]]
        vetor = [A, B]
        
        solucao = SistemaLinear.resolver(matriz, vetor)
        A1, A2 = solucao
        
        decomposicao.append({
            'tipo': 'linear_geral',
            'coeficiente': A1,
            'fator': fatores[0],
            'forma': f"{A1:.4g} / ({a1:.4g}x + {b1:.4g})"
        })
        decomposicao.append({
            'tipo': 'linear_geral',
            'coeficiente': A2,
            'fator': fatores[1],
            'forma': f"{A2:.4g} / ({a2:.4g}x + {b2:.4g})"
        })
    
    elif tipo_fatoracao == 'linear':
        # (Ax + B) / [(x - x₁)(x - x₂)] = A₁/(x - x₁) + A₂/(x - x₂)
        x1, x2 = info_fatoracao['raizes']
        
        # Sistema: A₁ + A₂ = A
        #         -A₁x₂ - A₂x₁ = B
        matriz = [[1, 1], [-x2, -x1]]
        vetor = [A, B]
        
        solucao = SistemaLinear.resolver(matriz, vetor)
        A1, A2 = solucao
        
        decomposicao.append({
            'tipo': 'linear',
            'coeficiente': A1,
            'raiz': x1,
            'forma': f"{A1:.4g} / (x - {x1:.4g})"
        })
        decomposicao.append({
            'tipo': 'linear',
            'coeficiente': A2,
            'raiz': x2,
            'forma': f"{A2:.4g} / (x - {x2:.4g})"
        })
    
    elif tipo_fatoracao == 'quadratico_complexo':
        # Mantém na forma (Ax + B) / (ax² + bx + c)
        a, b, c = info_fatoracao['coeficientes']
        decomposicao.append({
            'tipo': 'quadratico',
            'numerador': (A, B),
            'denominador': (a, b, c),
            'forma': f"({A:.4g}x + {B:.4g}) / ({a:.4g}x² + {b:.4g}x + {c:.4g})"
        })
    
    elif tipo_fatoracao == 'misto':
        # (Ax + B) / [(x² + p₁)(x² + p₂)] = (C₁x + D₁)/(x² + p₁) + (C₂x + D₂)/(x² + p₂)
        p1, p2 = info_fatoracao['valores']
        
        # Sistema 4x4:
        # C₁ + C₂ = 0
        # D₁ + D₂ = A
        # p₂C₁ + p₁C₂ = 0
        # p₂D₁ + p₁D₂ = B
        
        matriz = [
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [p2, 0, p1, 0],
            [0, p2, 0, p1]
        ]
        vetor = [0, A, 0, B]
        
        solucao = SistemaLinear.resolver(matriz, vetor)
        C1, D1, C2, D2 = solucao
        
        decomposicao.append({
            'tipo': 'quadratico',
            'numerador': (C1, D1),
            'denominador': (1, 0, p1),
            'forma': f"({C1:.4g}x + {D1:.4g}) / (x² + {p1:.4g})"
        })
        decomposicao.append({
            'tipo': 'quadratico',
            'numerador': (C2, D2),
            'denominador': (1, 0, p2),
            'forma': f"({C2:.4g}x + {D2:.4g}) / (x² + {p2:.4g})"
        })
    
    return decomposicao


def integraCadaTermo(termo):
    """
    Integra um termo da decomposição.
    
    Retorna:
        string com a integral do termo
    """
    if termo['tipo'] == 'linear':
        # ∫ A/(x - x₀) dx = A·ln|x - x₀|
        A = termo['coeficiente']
        x0 = termo['raiz']
        
        if abs(x0) < 1e-10:
            return f"{A:.4g}·ln|x|"
        else:
            return f"{A:.4g}·ln|x - ({x0:.4g})|"
    
    elif termo['tipo'] == 'linear_geral':
        # ∫ A/(ax + b) dx = (A/a)·ln|ax + b|
        A = termo['coeficiente']
        b, a = termo['fator'][0], termo['fator'][1]
        
        if abs(a) < 1e-10:
            return "erro: divisão por zero"
        
        coef_integral = A / a
        
        if abs(b) < 1e-10:
            return f"{coef_integral:.4g}·ln|{a:.4g}x|"
        else:
            return f"{coef_integral:.4g}·ln|{a:.4g}x + {b:.4g}|"
    
    elif termo['tipo'] == 'quadratico':
        # ∫ (Cx + D)/(ax² + bx + c) dx
        C, D = termo['numerador']
        a, b, c = termo['denominador']
        
        partes = []
        
        # Parte 1: termo com derivada (se C ≠ 0)
        if abs(C) > 1e-10:
            # ∫ C·x/(ax² + bx + c) dx relacionado com ∫ (2ax + b)/(ax² + bx + c) dx
            coef_ln = C / (2*a)
            if abs(b) < 1e-10:
                partes.append(f"{coef_ln:.4g}·ln|{a:.4g}x² + {c:.4g}|")
            else:
                partes.append(f"{coef_ln:.4g}·ln|{a:.4g}x² + {b:.4g}x + {c:.4g}|")
        
        # Parte 2: termo arctg (relacionado com a parte constante ajustada)
        delta = b**2 - 4*a*c
        
        if delta < 0:  # Raízes complexas
            # Completa quadrado: a(x + b/(2a))² + (c - b²/(4a))
            h = b / (2*a)
            k = (4*a*c - b**2) / (4*a)
            
            # Ajuste no numerador constante
            D_ajustado = D - C * b / (2*a)
            
            if abs(D_ajustado) > 1e-10 and abs(k) > 1e-10:
                coef_arctan = D_ajustado / (math.sqrt(a * k))
                arg = f"(x + {h:.4g}) / {math.sqrt(k/a):.4g}"
                partes.append(f"{coef_arctan:.4g}·arctan({arg})")
        
        return " + ".join(partes) if partes else "0"
    
    return "termo não reconhecido"


def calculaIntegral(numerador, denominador):
    """
    Função principal que calcula a integral completa.
    Coordena todas as etapas do processo.
    
    Retorna:
        dicionário com todos os passos e resultado
    """
    resultado = {
        'valido': False,
        'mensagem_validacao': '',
        'tipo_fatoracao': '',
        'info_fatoracao': {},
        'decomposicao': [],
        'integrais_parciais': [],
        'resultado_final': ''
    }
    
    # Passo 1: Verificar se a fração é válida
    valido, mensagem, den_expandido = verificaFracao(numerador, denominador)
    resultado['valido'] = valido
    resultado['mensagem_validacao'] = mensagem
    
    if not valido:
        return resultado
    
    # Passo 2: Identificar tipo de fatoração (usa denominador original se fatorado)
    tipo, info = identificaTipoFatoracao(denominador)
    resultado['tipo_fatoracao'] = tipo
    resultado['info_fatoracao'] = info
    
    # Passo 3: Decompor em frações parciais
    decomp = decompoeEmFracoesParciais(numerador, denominador, tipo, info)
    resultado['decomposicao'] = decomp
    
    # Passo 4: Integrar cada termo
    integrais = []
    for termo in decomp:
        integral = integraCadaTermo(termo)
        integrais.append(integral)
    resultado['integrais_parciais'] = integrais
    
    # Passo 5: Montar resultado final
    if integrais:
        resultado['resultado_final'] = " + ".join(integrais) + " + C"
    else:
        resultado['resultado_final'] = "C"
    
    return resultado


def imprimeResultado(resultado):
    """
    Imprime o resultado de forma formatada e organizada.
    """
    print("\n" + "="*70)
    print("RESULTADO DA INTEGRAÇÃO POR FRAÇÕES PARCIAIS")
    print("="*70)
    
    print(f"\n1. VALIDAÇÃO: {resultado['mensagem_validacao']}")
    
    if not resultado['valido']:
        return
    
    print(f"\n2. TIPO DE FATORAÇÃO: {resultado['tipo_fatoracao'].upper()}")
    print(f"   Descrição: {resultado['info_fatoracao']['descricao']}")
    if 'fatores_str' in resultado['info_fatoracao']:
        print(f"   Fatores: {resultado['info_fatoracao']['fatores_str']}")
    else:
        print(f"   Fatores: {resultado['info_fatoracao']['fatores']}")
    
    print(f"\n3. DECOMPOSIÇÃO EM FRAÇÕES PARCIAIS:")
    for i, termo in enumerate(resultado['decomposicao'], 1):
        print(f"   Termo {i}: {termo['forma']}")
    
    print(f"\n4. INTEGRAÇÃO DE CADA TERMO:")
    for i, integral in enumerate(resultado['integrais_parciais'], 1):
        print(f"   ∫ Termo {i} dx = {integral}")
    
    print(f"\n5. RESULTADO FINAL:")
    print(f"   ∫ f(x) dx = {resultado['resultado_final']}")
    print("\n" + "="*70 + "\n")


# ==============================================================================
# INTERFACE INTERATIVA
# ==============================================================================

def menuInterativo():
    """
    Interface interativa para o usuário inserir os valores.
    """
    print("\n" + "="*70)
    print("SISTEMA DE INTEGRAÇÃO POR FRAÇÕES PARCIAIS")
    print("="*70)
    print("\nFormato da integral: ∫ (numerador) / (denominador) dx")
    print("\nInsira os polinômios como expressões matemáticas:")
    print("Exemplos válidos:")
    print("  - 3x + 5")
    print("  - x^2 - 3x + 2")
    print("  - 10x^2 - 20")
    print("  - (x-1)(x+2)  [forma fatorada]")
    print("  - (3x-2)(x+5)  [forma fatorada]")
    print("="*70)
    
    while True:
        try:
            print("\n--- NUMERADOR ---")
            num_str = input("Digite o numerador (ex: 3x + 5): ")
            numerador = parsePolinomio(num_str)
            
            # Mostra interpretação
            if isinstance(numerador, dict) and numerador.get('fatorado'):
                print("Interpretado como: forma fatorada")
            else:
                print(f"Interpretado como: {Polinomio(numerador)}")
            
            print("\n--- DENOMINADOR ---")
            den_str = input("Digite o denominador (ex: x^2 - 3x + 2 ou (x-1)(x-2)): ")
            denominador = parsePolinomio(den_str)
            
            # Mostra interpretação
            if isinstance(denominador, dict) and denominador.get('fatorado'):
                print("Interpretado como: forma fatorada")
            else:
                print(f"Interpretado como: {Polinomio(denominador)}")
            
            # Calcular integral
            resultado = calculaIntegral(numerador, denominador)
            
            # Imprimir resultado
            imprimeResultado(resultado)
            
            # Perguntar se deseja continuar
            continuar = input("\nDeseja resolver outra integral? (s/n): ")
            if continuar.lower() != 's':
                print("\nEncerrando o programa. Até logo!")
                break
                
        except Exception as e:
            print(f"\nErro ao processar entrada: {e}")
            print("Tente novamente com o formato correto.")
            import traceback
            traceback.print_exc()


# ==============================================================================
# EXEMPLOS AUTOMÁTICOS
# ==============================================================================

def executarExemplos():
    """
    Executa exemplos automáticos para demonstração.
    """
    exemplos = [
        {
            'nome': 'a) ∫ (x-1)/(x²+2) dx',
            'numerador': 'x-1',
            'denominador': 'x^2+2'
        },
        {
            'nome': 'b) ∫ 1/(x²+2x+10) dx',
            'numerador': '1',
            'denominador': 'x^2+2x+10'
        },
        {
            'nome': 'c) ∫ (x⁴+9x-1)/(x+1)(x²+x+1) dx',
            'numerador': 'x^4+9x-1',
            'denominador': '(x+1)(x^2+x+1)'
        },
        {
            'nome': 'e) ∫ 2x/(x²+1) dx',
            'numerador': '2x',
            'denominador': 'x^2+1'
        },
        {
            'nome': 'f) ∫ 1/((x+2)(x-3)²) dx',
            'numerador': '1',
            'denominador': '(x+2)(x-3)(x-3)'
        },
        {
            'nome': 'g) ∫ 1/(x⁴-x²) dx',
            'numerador': '1',
            'denominador': 'x^4-x^2'
        },
        {
            'nome': 'h) ∫ (x+1)/(x(x+1)) dx',
            'numerador': 'x+1',
            'denominador': 'x(x+1)'
        },
        {
            'nome': 'Exemplo 1: Fatores lineares',
            'numerador': '3x + 5',
            'denominador': '(x-1)(x-2)'
        },
        {
            'nome': 'Exemplo 2: Fatorado geral',
            'numerador': '3x + 27',
            'denominador': '(3x-2)(x+5)'
        }
    ]
    
    for exemplo in exemplos:
        print(f"\n\n{'#'*70}")
        print(f"EXEMPLO: {exemplo['nome']}")
        print(f"Numerador: {exemplo['numerador']}")
        print(f"Denominador: {exemplo['denominador']}")
        print(f"{'#'*70}")
        
        try:
            num = parsePolinomio(exemplo['numerador'])
            den = parsePolinomio(exemplo['denominador'])
            
            resultado = calculaIntegral(num, den)
            imprimeResultado(resultado)
        except Exception as e:
            print(f"Erro ao processar este exemplo: {e}")
            import traceback
            traceback.print_exc()
        
        input("Pressione ENTER para continuar...")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    """
    Função principal - escolhe entre modo interativo ou exemplos.
    """
    print("\n" + "="*70)
    print("BEM-VINDO AO SISTEMA DE INTEGRAÇÃO POR FRAÇÕES PARCIAIS")
    print("="*70)
    print("\nEscolha uma opção:")
    print("1. Modo interativo (inserir valores manualmente)")
    print("2. Executar exemplos automáticos")
    print("3. Sair")
    
    opcao = input("\nOpção: ")
    
    if opcao == '1':
        menuInterativo()
    elif opcao == '2':
        executarExemplos()
    else:
        print("\nEncerrando o programa. Até logo!")


if __name__ == "__main__":
    main()