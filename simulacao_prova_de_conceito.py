"""
Prova de conceito do controle de proa - TCC CApA-Of 2026 (versao final)
1T Arthur dos Santos Sousa

Reproduz a simulacao da secao 5.9 (Apendice C) e gera os resultados da
secao 5.9.3:
  (a) otimizacao da razao de projeto r da regra de sintonia de
      Swider, Trybus e Stec (2023, Tabela 1);
  (b) sensibilidade a degradacao do atuador (limite de taxa e atraso).

Modelo: Nomoto de primeira ordem, K = 0,7 (graus/s)/grau e T = 60 s,
exemplo publicado por Swider, Trybus e Stec (2023). Nao representa o
navio estudado.

Requisitos: Python 3.12+, numpy e matplotlib (ver requirements.txt).
Uso: python simulacao_prova_de_conceito.py
Saidas (pasta ./saida, criada automaticamente; outra pasta pode ser
indicada pela variavel de ambiente TCC_SAIDA):
  fig09_proa_erro.png            Figura 11 do TCC
  fig10_ordem_perturbacao.png    Figura 12 do TCC
  fig11_compromisso_r.png        Figura 13 do TCC
  fig12_degradacao_atuador.png   Figura 14 do TCC
  resultados.txt                 numeros das Tabelas 14, 15 e 16 do TCC
"""
import os
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

# ---------------------------------------------------------------------------
# Parametros do cenario (identicos aos da secao 5.9 e do Apendice C)
# ---------------------------------------------------------------------------
K, T = 0.7, 60.0              # ganho [(graus/s)/grau] e constante de tempo [s]
DT, T_FIM = 0.05, 1000.0      # passo de integracao e duracao [s]
T_DEGRAU, REF = 30.0, 10.0    # instante e amplitude da mudanca de proa [s, graus]
PERT_INI, PERT_FIM, PERT_D = 480.0, 600.0, 0.10   # perturbacao [s, s, graus/s]
FAIXA = 0.5                   # faixa de assentamento [graus]
LIMITE_LEME = 35.0            # saturacao mecanica do leme [graus]

t = np.arange(0.0, T_FIM + DT / 2, DT)
ref = np.where(t >= T_DEGRAU, REF, 0.0)
pert = np.where((t >= PERT_INI) & (t < PERT_FIM), PERT_D, 0.0)

BASE = Path(__file__).resolve().parent
FIG = Path(os.environ.get("TCC_SAIDA", BASE / "saida"))
FIG.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Regras de sintonia
# ---------------------------------------------------------------------------
def sintonia_padrao(r=2.0, xi=1.0):
    """Regra-padrao de Fossen (2011), na forma da eq. (20) de
    Swider, Trybus e Stec (2023). Filtro de referencia Tf = Kd/Kp."""
    kp = r**2 / (K * T * xi**2)
    kd = (2 * r - 1) / K
    ki = r**3 / (10 * K * T**2 * xi**3)
    return kp, ki, kd, kd / kp


def sintonia_tabela1(r=2.0):
    """Regra da Tabela 1 (forma serie) de Swider, Trybus e Stec (2023).
    r = T/Tcl e a razao de projeto; filtro de referencia Tf = 2T/r."""
    kp = r * (r + 2) / (K * T)
    ki = r**2 / (K * T**2)
    kd = 2 * r / K
    return kp, ki, kd, 2 * T / r


# ---------------------------------------------------------------------------
# Simulacao (Euler explicito, condicoes iniciais nulas)
# ---------------------------------------------------------------------------
def simular(ganhos, taxa=None, atraso=0.0, t_servo=1.0):
    """Retorna proa, ordem de leme (delta_c) e posicao do leme (delta).

    taxa=None representa atuador ideal (delta = delta_c). Caso contrario,
    o atuador e um servo de primeira ordem (t_servo) com limite de taxa
    [graus/s], atraso puro [s] e saturacao em +/-35 graus.
    """
    kp, ki, kd, tf = ganhos
    n = len(t)
    psi_v, dc_v, d_v = np.zeros(n), np.zeros(n), np.zeros(n)
    psi = r_g = psi_f = integ = delta = 0.0
    n_atraso = int(round(atraso / DT))
    fila = [0.0] * (n_atraso + 1)
    for k in range(n):
        dpsi_f = (ref[k] - psi_f) / tf          # referencia filtrada
        e = psi_f - psi                         # erro de controle
        de = dpsi_f - r_g                       # derivada do erro
        dc = kp * e + ki * integ + kd * de      # ordem de leme
        dc = max(-LIMITE_LEME, min(LIMITE_LEME, dc))
        if taxa is None:
            delta = dc
        else:
            fila.append(dc)
            dc_atrasada = fila.pop(0)
            v = (dc_atrasada - delta) / t_servo
            v = max(-taxa, min(taxa, v))
            delta = max(-LIMITE_LEME, min(LIMITE_LEME, delta + v * DT))
        psi_v[k], dc_v[k], d_v[k] = psi, dc, delta
        integ += e * DT
        dr = (K * delta + pert[k] - r_g) / T   # Nomoto: T r' + r = K delta + d
        psi += r_g * DT
        r_g += dr * DT
        psi_f += dpsi_f * DT
    return psi_v, dc_v, d_v


# ---------------------------------------------------------------------------
# Indicadores
# ---------------------------------------------------------------------------
def indicadores(psi, dc, delta, limiar_div=0.2):
    erro = ref - psi                                   # erro de proa
    iae = float(np.sum(np.abs(erro)) * DT)
    janela = (t >= T_DEGRAU) & (t < PERT_INI)
    sobre = max(0.0, float(psi[janela].max() - REF))
    dentro = np.abs(psi - REF) <= FAIXA
    idx = np.where(janela)[0]
    t_ass = next((t[i] - T_DEGRAU for i in idx if dentro[i:idx[-1] + 1].all()), np.nan)
    jp = (t >= PERT_INI) & (t < PERT_FIM)
    desvio = float(np.abs(psi[jp] - REF).max())
    apos = np.where(t >= PERT_FIM)[0]
    recup = next((t[i] - PERT_FIM for i in apos if dentro[i:].all()), np.nan)
    atividade = float(np.sum(np.abs(np.diff(delta))))
    # tempo para o leme atingir 50% do pico da primeira ordem
    k0 = np.searchsorted(t, T_DEGRAU)
    kpico = k0 + int(np.argmax(dc[k0:k0 + int(60 / DT)]))
    k50 = k0 + int(np.argmax(delta[k0:] >= 0.5 * dc[kpico]))
    t50 = float(t[k50] - T_DEGRAU)
    dur_div = float(np.sum(np.abs(dc - delta) > limiar_div) * DT)
    return dict(iae=iae, assent=t_ass, sobre=sobre, desvio=desvio,
                recup=recup, atividade=atividade, t50=t50, dur_div=dur_div,
                pico=float(np.abs(dc).max()))


# ---------------------------------------------------------------------------
# Estilo das figuras (tinta neutra, grade discreta, duas series no maximo)
# ---------------------------------------------------------------------------
AZUL, LARANJA, VERDE = "#2a78d6", "#eb6834", "#1baf7a"
TINTA, TINTA2, GRADE = "#0b0b0b", "#52514e", "#e4e2dd"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif",
                                           "DejaVu Serif"],
    "font.size": 10, "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA,
    "xtick.color": TINTA2, "ytick.color": TINTA2, "axes.grid": True,
    "grid.color": GRADE, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 1.8,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})


class FormatadorVirgula(mticker.ScalarFormatter):
    """Rotulos de eixo com virgula decimal (padrao brasileiro)."""

    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace(".", ",")


def virgula(*eixos):
    for ax in eixos:
        ax.xaxis.set_major_formatter(FormatadorVirgula())
        ax.yaxis.set_major_formatter(FormatadorVirgula())


def sombrear_perturbacao(ax):
    ax.axvspan(PERT_INI, PERT_FIM, color="#f0efec", zorder=0, lw=0)


def main():
    linhas = []

    # ---------------- 5.9.1/5.9.2: comparacao das duas sintonias ----------
    g_pad, g_t1 = sintonia_padrao(2.0), sintonia_tabela1(2.0)
    psi_p, dc_p, d_p = simular(g_pad)
    psi_a, dc_a, d_a = simular(g_t1)
    m_p, m_a = indicadores(psi_p, dc_p, d_p), indicadores(psi_a, dc_a, d_a)
    linhas.append("== Comparacao das sintonias, r = 2 (Tabela 14 do TCC) ==")
    linhas.append("padrao = regra-padrao de Fossen (eq. 20 de Swider et al.); "
                  "tabela1 = regra da Tabela 1 de Swider et al. (2023)")
    for nome, g in (("padrao", g_pad), ("tabela1", g_t1)):
        linhas.append(f"ganhos {nome}: Kp={g[0]:.6f} Ki={g[1]:.8f} Kd={g[2]:.6f} Tf={g[3]:.1f}")
    for chave in ("iae", "assent", "sobre", "desvio", "recup", "atividade"):
        linhas.append(f"{chave:10s} padrao={m_p[chave]:9.3f}  tabela1={m_a[chave]:9.3f}")
    linhas.append("(assentamento contado a partir da mudanca de referencia, t = 30 s)")
    linhas.append(f"variacao IAE = {100 * (m_a['iae'] / m_p['iae'] - 1):.1f}%")
    linhas.append(f"variacao desvio = {100 * (m_a['desvio'] / m_p['desvio'] - 1):.1f}%")
    linhas.append(f"variacao atividade = {100 * (m_a['atividade'] / m_p['atividade'] - 1):.1f}%")
    linhas.append(f"acrescimo assentamento = {m_a['assent'] - m_p['assent']:.2f} s")

    # Figura: proa e erro
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.3, 4.6), sharex=True)
    for ax in (a1, a2):
        sombrear_perturbacao(ax)
    a1.plot(t, ref, color=TINTA2, lw=1.0, ls=(0, (4, 3)), label="Referência")
    a1.plot(t, psi_p, color=AZUL, label="Sintonia-padrão (Fossen)")
    a1.plot(t, psi_a, color=LARANJA, ls=(0, (6, 2)), label="Sintonia de Świder et al. (r = 2)")
    a1.set_ylabel("Proa (graus)")
    a1.legend(loc="lower right", fontsize=8.5)
    a2.plot(t, ref - psi_p, color=AZUL)
    a2.plot(t, ref - psi_a, color=LARANJA, ls=(0, (6, 2)))
    a2.axhline(0, color=TINTA2, lw=0.6)
    a2.set_ylabel("Erro de proa (graus)")
    a2.set_xlabel("Tempo (s)")
    a2.text(540, a2.get_ylim()[1] * 0.9, "perturbação", ha="center",
            va="top", fontsize=8, color=TINTA2)
    fig.align_ylabels()
    virgula(a1, a2)
    fig.savefig(FIG / "fig09_proa_erro.png")
    plt.close(fig)

    # Figura: ordem de leme e perturbacao
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.3, 4.2), sharex=True,
                                 gridspec_kw={"height_ratios": [2, 1]})
    for ax in (a1, a2):
        sombrear_perturbacao(ax)
    a1.plot(t, dc_p, color=AZUL, label="Sintonia-padrão (Fossen)")
    a1.plot(t, dc_a, color=LARANJA, ls=(0, (6, 2)), label="Sintonia de Świder et al. (r = 2)")
    a1.axhline(0, color=TINTA2, lw=0.6)
    a1.set_ylabel("Ordem de leme (graus)")
    a1.legend(loc="upper right", fontsize=8.5)
    a2.plot(t, pert, color=TINTA2, lw=1.4)
    a2.set_ylabel("Perturbação\n(graus/s)")
    a2.set_xlabel("Tempo (s)")
    fig.align_ylabels()
    virgula(a1, a2)
    fig.savefig(FIG / "fig10_ordem_perturbacao.png")
    plt.close(fig)

    # ---------------- 5.9.3 (a): otimizacao da razao de projeto ----------
    grade_r = np.round(np.arange(1.0, 4.0001, 0.05), 2)
    res = []
    for r in grade_r:
        psi, dc, d = simular(sintonia_tabela1(r))
        m = indicadores(psi, dc, d)
        res.append((r, m["iae"], m["atividade"], m["desvio"], m["sobre"], m["pico"]))
    res = np.array(res)
    i2 = int(np.argmin(np.abs(res[:, 0] - 2.0)))
    iae2, at2 = res[i2, 1], res[i2, 2]
    linhas.append("\n== Varredura de r na regra da Tabela 1 de Swider et al. (Tabela 15 e Figura 13 do TCC) ==")
    linhas.append("r     IAE      atividade  desvio  sobressinal  pico_ordem")
    for row in res:
        if abs(row[0] * 4 - round(row[0] * 4)) < 1e-9:      # passos de 0,25
            linhas.append(f"{row[0]:4.2f}  {row[1]:7.1f}  {row[2]:7.2f}   {row[3]:5.2f}   "
                          f"{row[4]:6.3f}      {row[5]:5.2f}")
    otimos = {}
    for lam in (0.5, 1.0, 2.0):
        J = res[:, 1] / iae2 + lam * res[:, 2] / at2
        j = int(np.argmin(J))
        otimos[lam] = res[j]
        linhas.append(f"lambda={lam:3.1f}: r*={res[j, 0]:.2f}  IAE={res[j, 1]:.1f}  "
                      f"atividade={res[j, 2]:.2f}  J={J[j]:.3f}")

    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    ax.plot(res[:, 2], res[:, 1], color=AZUL, lw=1.6, zorder=2)
    for rr in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        j = int(np.argmin(np.abs(res[:, 0] - rr)))
        ax.plot(res[j, 2], res[j, 1], "o", ms=5, color=AZUL, zorder=3)
        ax.annotate(f"r = {str(rr).replace('.', ',')}", (res[j, 2], res[j, 1]),
                    textcoords="offset points", xytext=(7, 4), fontsize=8.5, color=TINTA2)
    marcadores = {0.5: "s", 1.0: "D", 2.0: "^"}
    for lam, row in otimos.items():
        ax.plot(row[2], row[1], marcadores[lam], ms=9, mfc="white", mec=LARANJA,
                mew=1.8, zorder=4,
                label=f"ótimo para λ = {str(lam).replace('.', ',')} "
                      f"(r = {row[0]:.2f})".replace(".", ","))
    ax.plot(m_p["atividade"], m_p["iae"], "X", ms=8, color=TINTA2, zorder=4,
            label="sintonia-padrão (Fossen)")
    ax.set_xlabel("Atividade do leme Aδ (graus)")
    ax.set_ylabel("IAE (grau·s)")
    ax.legend(loc="upper right", fontsize=8.5)
    virgula(ax)
    fig.savefig(FIG / "fig11_compromisso_r.png")
    plt.close(fig)

    # ---------------- 5.9.3 (b): degradacao do atuador ---------------------
    cenarios = {
        "nominal": (2.3, 0.5),         # 65 graus em 28 s (SOLAS) e 0,5 s de atraso
        "degradado": (0.5, 3.0),
        "muito degradado": (0.2, 8.0),
    }
    linhas.append("\n== Degradacao do atuador, regra da Tabela 1 de Swider et al., r = 2 (Tabela 16 e Figura 14 do TCC) ==")
    linhas.append("cenario            IAE     desvio  t50[s]  t_div>0,2[s]  atividade")
    series = {}
    for nome, (taxa, atraso) in cenarios.items():
        psi, dc, d = simular(g_t1, taxa=taxa, atraso=atraso)
        m = indicadores(psi, dc, d)
        series[nome] = (psi, dc, d)
        linhas.append(f"{nome:17s} {m['iae']:7.1f}  {m['desvio']:5.2f}  {m['t50']:6.1f}  "
                      f"{m['dur_div']:8.1f}     {m['atividade']:6.2f}")

    psi_n, dc_n, d_n = series["nominal"]
    psi_d, dc_d, d_d = series["muito degradado"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.3, 4.8))
    a1.plot(t, psi_n, color=AZUL, label="Atuador nominal")
    a1.plot(t, psi_d, color=LARANJA, ls=(0, (6, 2)), label="Atuador muito degradado")
    a1.plot(t, ref, color=TINTA2, lw=1.0, ls=(0, (4, 3)), label="Referência")
    a1.set_xlim(0, 300)
    a1.set_ylabel("Proa (graus)")
    a1.set_xlabel("Tempo (s)")
    a1.legend(loc="lower right", fontsize=8.5)
    a1.set_title("(a) Proa", loc="left", fontsize=10)
    jan = (t >= 25) & (t <= 70)
    a2.plot(t[jan], dc_d[jan], color=TINTA2, lw=1.0, ls=(0, (4, 3)),
            label="Ordem de leme (muito degradado)")
    a2.plot(t[jan], d_n[jan], color=AZUL, label="Posição do leme (nominal)")
    a2.plot(t[jan], d_d[jan], color=LARANJA, ls=(0, (6, 2)),
            label="Posição do leme (muito degradado)")
    a2.set_ylabel("Ângulo de leme (graus)")
    a2.set_xlabel("Tempo (s)")
    a2.set_ylim(-0.3, 1.75)
    a2.legend(loc="upper right", fontsize=8.5)
    a2.set_title("(b) Ordem e posição do leme", loc="left", fontsize=10)
    virgula(a1, a2)
    fig.tight_layout(h_pad=1.2)
    fig.savefig(FIG / "fig12_degradacao_atuador.png")
    plt.close(fig)

    texto = "\n".join(linhas)
    (FIG / "resultados.txt").write_text(texto + "\n", encoding="utf-8")
    print(texto)


if __name__ == "__main__":
    main()
