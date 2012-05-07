import pylab
import parameters
import matplotlib

def plot(filesuffix='.png'):

    matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

    scale = 1e6
    d = parameters.delta * scale
    r = parameters.fieldWidth * scale
    h = parameters.featureDepth * scale
    w = parameters.trenchWidth * scale

    fig = pylab.figure()
    ax = fig.add_subplot(111, aspect='equal')

    extra = 0

    ax.add_patch(pylab.Rectangle((-40, -h), 200, h + extra, alpha=0.3, color='black', ec='none'))
    ax.add_patch(pylab.Rectangle((-40, extra), 200, d - extra, alpha=0.2, color='black', ec='none'))
    ax.add_patch(pylab.Rectangle((-40, d), 200, 20, alpha=0.1, color='black', ec='none'))
    ax.add_patch(pylab.Rectangle((-40, d + 20), 200, 20, alpha=0.05, color='black', ec='none'))
    ax.add_patch(pylab.Rectangle((-40, -h - 40), 200, 40, alpha=0.5, color='black', ec='none'))

    pylab.plot((-r / 2, -r / 2, -w / 2, -w / 2, w / 2, w / 2, r / 2, r / 2, -r / 2),
               (d, 0, 0, -h, -h, 0, 0, d, d), 'k', lw=3, alpha=0.5)

    pylab.ion()
    pylab.xticks((-h / 2, h / 2), (r'$-a_S^{} / P$', r'$a_S^{} / P$'))
    pylab.yticks((-56, 0, 150, 170), (r'$-h$', 0, r'$\delta$', r'$L$'))

    # ax2 = pylab.twinx(aspect='equal')
    # pylab.yticks((-56, 0, 150), (r'$-h$', r'0',r'$\delta$'))
    # pylab.ylim(ymin=-0.00007 * scale)
    # pylab.ylim(ymax=0.000180 * scale)
    # pylab.xticks((-0.00004 * scale, 0, 0.00004 * scale))
    # pylab.xlim(xmax=0.00008 * scale)
    # pylab.xlim(xmin=-0.00004 * scale)

    pylab.axes(ax, aspect='equal')
    pylab.ylim(ymin=-0.00007 * scale)
    pylab.ylim(ymax=0.000180 * scale)
    pylab.xlim(xmax=0.00008 * scale)
    pylab.xlim(xmin=-0.00004 * scale)
    pylab.xlabel(r'$x$')
    pylab.ylabel(r'$z$')

    pylab.text(-25, 70, 
                r'\['
                r'\begin{split}'
                r'    \partial_z^2 \eta_{\text{Eff}}^{} &= 0 \\'
                r'    \partial_t C_{\text{Supp}}^{} &= D_{\text{Supp}}^{} \partial_z^2 C_{\text{Supp}}^{} \\'
                r'    \partial_t C_{\text{Cu}}^{} &= D_{\text{Cu}}^{} \partial_z^2 C_{\text{Cu}}^{}'
                r'\end{split}'
                r'\]', fontsize=8)

    pylab.text(-24, -16,
                r'\['
                r'\begin{split}'
                r'    \partial_t \theta =& k^+ C_{\text{Supp}}^{} \left(1 - \theta\right) \\'
                r'    &- k^- \theta v'
                r'\end{split}'
                r'\]', rotation='vertical', fontsize=8)
    pylab.text(-2, -25, r'Feature', rotation='vertical', fontsize=8)

    pylab.text(5, -50,
                r'\['
                r'\begin{split}'
                r'    \partial_z^2  \eta_{\text{Eff}}^{} =& \frac{\rho n F}{\Omega} v \Theta \\'
                r'    \partial_t C_{\text{Supp}}^{} =& \partial_z^2 C_{\text{Supp}}^{} \\'
                r'    - \Gamma k^+ & C_{\text{Supp}}^{} \left(1-\theta\right)  \Theta\\'
                r'    \partial_t C_{\text{Cu}}^{} =& D_{\text{Cu}}^{} \partial_z^2 C_{\text{Cu}}^{}  - \frac{v}{\Omega} \Theta '
                r'\end{split}'
                r'\]', rotation='horizontal', fontsize=8)
    pylab.text(-10, 105, r'Electrolyte', fontsize=8)

    pylab.text(-27, 152, 
                r'$C_{\text{Supp}}^{} = C_{\text{Supp}}^{\infty}, C_{\text{Cu}}^{} = C_{\text{Cu}}^{\infty}$', fontsize=8)
    pylab.text(-10, 162, r'Boundary Layer', fontsize=8)

    pylab.text(-27, 172, 
                r'$\eta_{\text{Eff}}^{} = E_{\text{Applied}}^{}$', fontsize=8)
    pylab.text(15, 173, r'Reference Electrode', fontsize=8)

    pylab.text(-30, -68, r'$\Theta \left(z\right)=\frac{P}{a_F^{}} + \delta \left(z\right) \left(\frac{a_S^{}}{a_F^{}} - 1\right) + \delta\left(z+h\right)$', fontsize=8)
    ##pylab.text(57, 87, r'$L$', fontsize=8)

    ##from matplotlib.patches import FancyArrowPatch
    ##ax.add_patch(FancyArrowPatch((55, -3),(55, 175),arrowstyle='<->', lw=2, mutation_scale=20))

    pylab.savefig('schematic' + filesuffix)

if __name__ == '__main__':
    plot()
    pylab.show()
