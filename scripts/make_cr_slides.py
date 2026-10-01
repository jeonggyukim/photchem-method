"""Slides: how radiation, cosmic rays, photochemistry and gas dynamics couple in the
TIGRESS++ framework, and the self-confinement picture of cosmic-ray transport.  The
explanations are in the presenter notes.

    python make_cr_slides.py [output.pptx]

Writes ~/Dropbox/Research-Presentation/2026-10-Pisa-Salt-of-SF/cr-coupling.pptx unless
an output path is given.  Figures are read from ../figures (git-ignored); make them
first with rad_cr_chem_coupling.py beside this script:

    python rad_cr_chem_coupling.py 3d2 talk --mhd --cr-left
    python rad_cr_chem_coupling.py 3d2 talk --mhd --cr-left --crpic=A
    python rad_cr_chem_coupling.py 3d2 talk --mhd --cr-left --crpic=B
    python rad_cr_chem_coupling.py 3d2 paper
"""
import os
import sys

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', 'figures')
OUT = (sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    '~/Dropbox/Research-Presentation/2026-10-Pisa-Salt-of-SF/cr-coupling.pptx'))

NAVY = RGBColor(0x1F, 0x3A, 0x5F)
GREY = RGBColor(0x40, 0x40, 0x40)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]
W, H = prs.slide_width, prs.slide_height


def title(slide, text):
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, W, Inches(0.9))
    bar.fill.solid()
    bar.fill.fore_color.rgb = NAVY
    bar.line.fill.background()
    tf = bar.text_frame
    tf.margin_left = Inches(0.5)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(28)
    p.font.bold = True
    p.font.color.rgb = WHITE


def figure(slide, name, top=Inches(1.0), max_h=Inches(6.4), max_w=Inches(12.9)):
    path = os.path.join(FIG, name)
    pic = slide.shapes.add_picture(path, 0, top, height=max_h)
    if pic.width > max_w:
        ratio = max_w / pic.width
        pic.width, pic.height = int(max_w), int(pic.height*ratio)
    pic.left = int((W - pic.width)/2)
    return pic


def notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text.strip()


def bullets(slide, items, top=Inches(1.15), size=18):
    tb = slide.shapes.add_textbox(Inches(0.6), top, W - Inches(1.2), H - top - Inches(0.3))
    tf = tb.text_frame
    tf.word_wrap = True
    for k, (level, s) in enumerate(items):
        p = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
        p.text = ('• ' if level == 0 else '– ') + s
        p.level = level
        p.font.size = Pt(size - 2*level)
        p.font.color.rgb = GREY
        p.space_after = Pt(6)
        p.alignment = PP_ALIGN.LEFT


# 1. the coupled system ---------------------------------------------------------------
s = prs.slides.add_slide(BLANK)
title(s, 'Radiation, cosmic rays, photochemistry and gas dynamics, coupled')
figure(s, 'rad_cr_chem_coupling_3d2_talk_mhd_crleft.png')
notes(s, """
Four modules, each with the arrows of what it hands the others.

Radiation transfer. Point sources use adaptive rays that split as they spread
(Abel & Wandelt 2002); the diffuse field uses parallel rays in fixed directions, with
the external background and scattering. Both cross the same clumpy cloud here, and each
ray fades with the optical depth it has crossed. From the specific intensity we take
the energy density E_rad = (1/c) ∫∮ I_ν dΩ dν and the flux F_rad = ∫∮ I_ν n̂ dΩ dν.

Radiation -> chemistry: band energy densities E_LyC, E_FUV and the shielded
(pseudo) energy density Ẽ_H2 drive photoionization, photodissociation and photoelectric
heating. These are 4π/c times the band mean intensities J of Kim et al. (2023).

Chemistry -> radiation: the opacity χ_ν (dust, HI, H2) and the self-shielding
columns N_H2, N_C, N_CO.

Cosmic-ray transport: the two-moment equations of Armillotta et al. (2021).
- v_m is the reduced speed of light (Jiang & Oh 2018): the maximum signal speed of the
  scheme, well below c but above every other speed. It sets the time step without
  changing the steady flux; when it is large enough, ∂F_c/∂t is negligible and the flux
  equation reduces to diffusion plus streaming. v_p is the real particle speed (≈ c).
- Q is the CR energy injected per unit volume and time by supernovae (10 per cent of
  E_SN in Armillotta et al. 2021). Λ_coll is the collisional loss (pion production,
  Coulomb and ionization).
- G ≡ σ_tot·[F_c − (4/3) v e_c] is the momentum exchanged between CRs and gas, the
  force density CRs exert on the gas. The bracket is the CR flux in the gas frame
  ((4/3) e_c v = (e_c + P_c) v is the enthalpy carried with the gas). -G appears in the
  CR flux equation and +G in the gas momentum equation. In steady state G = −∇P_c.
- In the energy equation −(v + v_s)·G splits into −v·G, work done by CR pressure on
  the moving gas, and −v_s·G, energy CRs put into the waves while streaming, which the
  damping turns into gas heat (streaming heating Γ_st = |v_s·∇P_c|). The gas total
  energy gains (v + v_s)·G; the momentum equation (+G) gives the kinetic energy v·G;
  so the internal energy gains only v_s·G. The work term heats the gas only later,
  through compression or shocks of the flow it drives. Hence the CR -> gas arrow
  carries only the force G, and Γ_st goes on the CR -> chemistry arrow as heating.
- σ_tot⁻¹ = σ_∥⁻¹ + v_A,i (P_c + e_c)/|B̂·∇P_c|: diffusion past the waves plus streaming
  at the ion Alfvén speed, written as an effective scattering. In steady state the flux
  equation gives F_∥ − (4/3) v_∥ e_c = −∇_∥P_c/σ_tot, so 1/σ acts as a diffusivity.
  Streaming carries F_st = (e_c + P_c) v_A,i down the gradient; writing it as
  |∇_∥P_c|/σ_st gives σ_st⁻¹ = v_A,i (P_c + e_c)/|B̂·∇P_c| = v_A,i L_c, with
  L_c = 4P_c/|∇_∥P_c| the gradient length. The diffusive and streaming fluxes add,
  so the inverse σ's add. κ_st depends on the gradient itself (nonlinear, → ∞ as the
  gradient flattens); the reduced speed of light v_m caps the flux.
- ∇P_c is a gradient of a scalar: the tensor form ∇·P_c reduces to ∇P_c for
  isotropic pressure P_c = e_c/3.

Radiation -> cosmic rays (dashed): inverse-Compton losses of CR electrons need E_rad
(Linzer et al. 2025). The power is P = (4/3) σ_T c γ² U_rad. At fixed energy σ_T ∝ m⁻²
and γ² ∝ m⁻², so P ∝ m⁻⁴: a proton loses (1836)⁴ ≈ 10¹³ times less than an electron. A
1 GeV electron in 1 eV cm⁻³ loses its energy in ≈ 3×10⁸ yr. TIGRESS transports protons
only, hence dashed.

Cosmic rays -> chemistry: ionization ξ_cr ∝ e_c, ionization heating Γ_cr, and
streaming heating Γ_st.
Chemistry -> cosmic rays: x_i, x_n and T set σ_∥ (through wave damping) and v_A,i.

Gas dynamics: star formation, gravity, galactic shear, stellar feedback. It takes the
radiation force f_rad, the CR force G and the net heating nΓ − n²Λ, and supplies ρ, v,
B to all the others. Implementation note: the Tigris cr_mg source step adds the energy
CRs lose to the gas total energy (u(IEN) −= Δe_c) and −ΔF_c/v_m to its momentum, so
streaming heating already reaches the gas in an MHD run; it is not yet in the
chemistry's nΓ, which matters in post-processing where the source coupling is off.

Photochemistry and thermodynamics: dx_i/dt = C_i − D_i x_i and de/dt = nΓ − n²Λ, with
rates depending on n, T, x_s, Z'_g, Z'_d, E and ξ_cr (Kim et al. 2023, Eqs. 12–14; Λ
also on |dv/dr| through CO trapping).
""")

# 2. resonant scattering ----------------------------------------------------------------
s = prs.slides.add_slide(BLANK)
title(s, 'Self-confinement: CRs scatter off the waves they drive')
figure(s, 'rad_cr_chem_coupling_3d2_talk_mhd_crleft_crA.png')
notes(s, """
The cartoon: a CR streams along B and meets two Alfvén-wave packets. The first is long
and matches the CR's step per gyration before scattering; the CR's pitch angle changes
there and its helix tightens. The second packet is shorter and matches the new, smaller
step: a CR at smaller pitch-angle cosine resonates with shorter waves, k = 1/(μ r_L).

Why resonance. The wave's δB is perpendicular to B, and the force is
F = (q/c) v × δB. The part that changes the pitch angle is v_⊥ × δB, along B. v_⊥
rotates at the gyrofrequency Ω/γ; moving along the field at v_∥ through a nearly static
wave (v_A ≪ c), the CR sees δB rotate at k v_∥. If the two differ, the phase between
v_⊥ and δB drifts and the force averages to zero. If k v_∥ = Ω/γ — one wavelength per
gyration, λ = 2π μ r_L — the phase stays fixed and the kicks add, like pushing a swing
at its natural frequency. Real waves come in packets with random phases, so each
resonant packet gives a coherent kick Δμ ~ δB/B in a random direction: a random walk in
pitch angle with rate ν ~ Ω (δB/B)², which is what sets σ_∥. Alfvén waves are circularly
polarized, so each resonates with particles moving one way along B. The streaming
instability is the same resonance run in reverse.

Scale. For a 1 GeV proton (pc = 1.7 GeV) in 3 μG, r_L = pc/(eB) ≈ 2×10¹² cm ≈ 0.1 AU.
The waves are ~10⁸ times smaller than a pc-scale cell, so σ_∥ comes from a local
growth–damping balance rather than resolved waves.

Resonant momentum. A wave of wavenumber k resonates with every CR with p μ = p₁ = eB/(ck)
= m_p Ω₀/k. Below p₁ the gyroradius is shorter than 1/k; the CR gyrates many times per
wavelength, sees an adiabatic change, and is not scattered. n₁ = 4π ∫_{p₁}^∞ p² F dp
counts the CRs that can drive and be scattered by those waves. Because the spectrum is
steep (F ∝ p^−4.7), n₁ is dominated by CRs just above p₁: each energy has its own wave
population, which is why transport is energy dependent (Armillotta et al. 2025) and why
self-confinement fails at very high energy. n₁ = 1.1×10⁻¹⁰ (e_c/eV) cm⁻³ for a fixed
spectral slope, so n₁/P_c is constant.

The 90° problem. Crossing μ = 0 needs k → ∞, where there is little wave power and
strong damping; quasilinear scattering stalls there. Mirroring, resonance broadening,
finite wave speed and transit-time damping of compressive modes carry CRs across. The
two-moment scheme does not track μ: one σ_∥ per energy group, averaged over pitch angle
and the resonant band (the √π/16 in A22 Eq. 16 comes from those integrals).
""")

# 3. the self-confinement loop -----------------------------------------------------------
s = prs.slides.add_slide(BLANK)
title(s, 'Self-confinement loop: growth, scattering, damping')
figure(s, 'rad_cr_chem_coupling_3d2_talk_mhd_crleft_crB.png')
notes(s, """
1. ∇P_c, the free energy. CRs are injected by supernovae, clustered near the midplane,
and lost in dense gas (Λ_coll) and at the boundaries, so P_c cannot be uniform. CRs
gyrate with r_L ~ 10¹² cm and move freely along B but hardly across it (σ_⊥ = 10 σ_∥ in
A22), so the excess drains along field lines and the gradient that drives transport is
B̂·∇P_c. Only it drives streaming and waves; the streaming velocity
v_s = −v_A,i B̂ sign(B̂·∇P_c) points along B, down the gradient.

2. Streaming, v_D > v_A,i -> waves δB: the streaming instability. In the wave frame CRs
drifting faster than v_A,i are anisotropic; resonant CRs hand momentum to forward Alfvén
waves, which grow at Γ_grow ∝ Ω₀ (n₁/n_i)(v_D/v_A,i − 1).

3. Waves -> scattering -> σ_∥. Pitch-angle scattering at ν ~ Ω(δB/B)²; σ_∥ ~ ν/v_p².
Stronger waves, more scattering, slower diffusion (κ_∥ ∝ 1/σ_∥).

4. The dashed arrow: scattering isotropizes CRs in the wave frame, removing the
anisotropy that drives the waves, so v_D → v_A,i and growth switches off. CRs stream at
≈ v_A,i instead of c.

5. Wave damping: ion–neutral and nonlinear Landau. Growth = damping fixes δB/B and so
σ_∥ (A22 Eqs. 16–17):
  σ_∥,NLL ∝ |B̂·∇P_c|^1/2 T^−1/4 n_i^−1/4   (hot ionized gas; thermal ions absorb the
                                            beat of two waves; rate ∝ T^1/2)
  σ_∥,IN  ∝ |B̂·∇P_c| n_i^−1/2 n_n^−1       (partly neutral gas; neutrals drag the ions
                                            that carry the wave; rate ∝ n_n)
The full forms carry |B̂·∇P_c| n₁ / (v_A,i P_c); n₁ ∝ P_c makes the 1/P_c drop out.

Neutrals: more neutrals damp the waves harder, δB is smaller, σ_∥ is smaller, κ_∥ ∝ n_n
is larger: faster diffusion. Mostly neutral gas also has a large v_A,i, so streaming is
faster too. Armillotta et al. (2021) find strong diffusion and nearly uniform P_c in most
of the ISM mass; fast diffusion flattens the gradient, which lowers σ_∥,IN further.

Ions: n_i^−1/2 comes from 1/n_i (growth: fewer ions to push, faster growth, more
scattering) times 1/v_A,i ∝ n_i^1/2 (drive: larger v_A,i, smaller drift excess). At fixed
n_H, σ_∥,IN ∝ x_i^−1/2 n_H^−3/2: lower ionization means slower diffusion but faster
streaming (v_A,i ∝ x_i^−1/2). The ion fraction from the chemistry enters three times: wave
growth, wave speed, and damping through n_n — hence x_i, x_n, T on the chemistry ->
CR arrow.

The damping heats the gas: the energy CRs lose to the waves is the streaming heating
v_s·G. σ_∥ enters σ_tot, which sets G.

One line: CRs streaming faster than the ion Alfvén speed excite the waves that scatter
them; scattering slows them toward v_A,i; ion–neutral and nonlinear Landau damping set
the wave amplitude, so the gas composition and temperature decide how fast CRs move.
""")

# 4. energy losses ------------------------------------------------------------------------
s = prs.slides.add_slide(BLANK)
title(s, 'Cosmic-ray energy losses')
bullets(s, [
    (0, 'Protons: ionization/Coulomb at low energy, pion production above ~1 GeV'),
    (1, 'Bethe: dE/dt ∝ 1/β, so t_loss ∝ E^3/2 (non-relativistic), ∝ E (relativistic)'),
    (1, 'n_H = 1 cm⁻³: t_ion ≈ 2×10⁵ yr at 10 MeV, ≈ 2×10⁸ yr at 1 GeV; '
        't_π ≈ 7×10⁷ yr'),
    (1, 'all ∝ 1/n_H: a 10 MeV proton lasts ~20 yr at n_H = 10⁴ cm⁻³'),
    (0, 'Electrons: inverse Compton and synchrotron'),
    (1, 'P_IC ∝ σ_T γ² U_rad ∝ m⁻⁴ at fixed E: 10¹³ times weaker for protons'),
    (1, '1 GeV electron, U_rad = 1 eV cm⁻³: t_IC ≈ 3×10⁸ yr'),
    (0, 'Ionization of dense gas comes from CRs below ~100 MeV'),
    (1, 'the ones that lose energy fastest, so ξ_cr falls with column inside clouds'),
    (1, 'one ~GeV group is transported; ξ_cr ∝ e_c with a fixed spectrum, attenuated by N_eff'),
])
notes(s, """
Ionization and Coulomb losses follow the Bethe stopping power, dE/dx ∝ Z²/β², so
dE/dt = v dE/dx ∝ 1/β: slow particles lose faster and have less to lose. Times above
use the Mannheim & Schlickeiser (1994) fit |dE/dt| ≈ 1.8×10⁻⁷ [2β²/(10⁻⁶ + 2β³)] n_H eV s⁻¹
(order of magnitude, factor ~2). Pion production has a threshold near 280 MeV and a
nearly constant cross-section. In ionized gas Coulomb losses on free electrons are a few
times the neutral-gas ionization losses.

Inverse Compton: P = (4/3) σ_T c γ² U_rad; σ_T ∝ (m_e/m)² and γ = E/mc² give P ∝ m⁻⁴ at
fixed E. For protons IC is irrelevant; their losses (Λ_coll) do not need U_rad.

Low-energy CRs: most of ξ_cr comes from CRs below ~100 MeV (Padovani et al. 2009,
2018), the ones with the shortest loss times, so ξ_cr drops with column inside dense
clouds. Tigris transports one ~1 GeV group with Λ_coll at 1 GeV; the chemistry takes
ξ_cr ∝ e_c with a fixed spectrum and applies a column attenuation (the -shld_cr path,
N_eff). Multi-group cr_mg transport (Armillotta et al. 2025) would let the low-energy
spectrum steepen on its own.
""")

# 5. the paper version -------------------------------------------------------------------
s = prs.slides.add_slide(BLANK)
title(s, 'Backup: post-processing version (paper figure)')
figure(s, 'rad_cr_chem_coupling_3d2_paper.png')
notes(s, """
The post-processing version: no gas dynamics box; radiation transfer, cosmic-ray
transport and photochemistry iterate on a fixed MHD snapshot. Labels name the processes
on each arrow. Made with rad_cr_chem_coupling.py 3d2 paper.
""")

prs.save(OUT)
print('wrote', OUT)
