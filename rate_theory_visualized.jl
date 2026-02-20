### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> image = "https://external-content.duckduckgo.com/iu/?u=https%3A%2F%2Fhoomd-blue.readthedocs.io%2Fen%2Fv5.0.1%2F_images%2Ftutorial_01-Introducing-Molecular-Dynamics_02-Initializing-a-Random-System_42_0.png&f=1&nofb=1&ipt=30e3e783029319fc4717b96e151aa316e94577a36ef4c0fe3395c11c3bf3b80a"
#> language = "en-US"
#> title = "Rate Theory Visualized"
#> date = "2026-02-15"
#> tags = ["Rate theory"]
#> description = "This document illustrates rate theory using interactivity and animations, as well as derivations, in order to confer a deep undersanding of rate theory."
#> 
#>     [[frontmatter.author]]
#>     name = "Nathaniel Thomas"
#>     url = "https://nathanielwthomas.xyz"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ f5bc5f70-0a80-11f1-a259-8be415f4d105
using PlutoUI, Handcalcs, CairoMakie, WGLMakie, Symbolics, Match, Latexify, Random, PhysicalConstants.CODATA2018, Unitful, GeometryBasics, Makie

# ╔═╡ 6d0c96f2-4dfd-43a1-86af-e6aaa9b04f0c
TableOfContents()

# ╔═╡ 4b5af7d2-3894-4220-8f2d-dec15268305d
md"""
# Rate Theory Visualized

This document provides a detailed overview of the mathematics behind rate theory.

This theory can be used to develop a theoretical framework for many phenomena, including corrosion, as well as void nucleation.
"""

# ╔═╡ 0dc3df3c-0999-4f9c-8e68-4b71a51d00ed
wglmakieswitch = @bind wglmakie Switch(default=true);

# ╔═╡ 45f0efff-bc20-49d3-8442-3948d82b6441
md"""
Use interactive visualizations: $wglmakieswitch
"""

# ╔═╡ 8b4c2f8f-6027-43f9-8bc2-d2fbd7cbf3ea
md"""
## Einstein relation

One of the easiest ways to represent diffusion is by the Einstein relation. This assumes that a particle diffusises randomly with no bias in any direction. It simply takes unit-length single steps in a random direction in 3D space. The graph below uses the following two sliders to render a given number of particles and steps to show what the Einstein relation looks like.
"""

# ╔═╡ 1be863ee-0514-4de1-a0c4-9460cf2852c9
step_slider = @bind steps PlutoUI.Slider(1:250; default=3, show_value=true);

# ╔═╡ a0a819cd-7d17-4d61-a88d-a2a9fdd3e8ab
particle_slider = @bind particles PlutoUI.Slider(1:100; default=1, show_value=true);

# ╔═╡ c8a756c6-14d3-4bb9-ab4f-e794a0653344
md"""
Adjust the sliders to choose the number of steps and particles.

Number of steps = $step_slider $(steps == 1 ? "step" : "steps")

Number of particles = $particle_slider $(particles == 1 ? "particle" : "particles").
"""

# ╔═╡ 48a4e303-921c-4b59-93a7-633f2f2cffb8
# This is just plotting code, no need to show.
begin
	if wglmakie
		WGLMakie.activate!()
	else 
		CairoMakie.activate!()
	end
	
	f = Figure()
	a = Axis3(f[1,1], title="Random Walk (Einstein Relation)")
	for particle in 1:particles
		x = [0]
		y = [0]
		z = [0]
		for step in 1:steps
			# There are six possible directions, choose one
			direction = rand(1:6)
			@match direction begin
				1 => begin
					append!(x, x[end] + 1)
					append!(y, y[end])
					append!(z, z[end])
				end
				2 => begin
					append!(x, x[end] - 1)
					append!(y, y[end])
					append!(z, z[end])
				end
				3 => begin
					append!(x, x[end])
					append!(y, y[end] + 1)
					append!(z, z[end])
				end
				4 => begin
					append!(x, x[end])
					append!(y, y[end] - 1)
					append!(z, z[end])
				end
				5 => begin
					append!(x, x[end])
					append!(y, y[end])
					append!(z, z[end] + 1)
				end
				6 => begin
					append!(x, x[end] + 1)
					append!(y, y[end])
					append!(z, z[end] - 1)
				end
			end
		end
		lines!(a, x, y, z)
	end
	f
end

# ╔═╡ 14b9e9d3-a71b-42a4-9942-94ebfe63ccb4
@variables R X Y Z;

# ╔═╡ db101921-c3dd-4d69-956a-5bc6650dd273
md"""
We can express the distance of a particle after $n$ steps from the origin using: 

$$R^2 = X^2 + Y^2 + Z^2 \tag{1}$$
"""

# ╔═╡ 3e1815e7-8306-4aa3-ada8-f113d2f4785f
md"""
where $X$, $Y$, and $Z$ are components along the $x$, $y$, and $z$ axis.

If you repeat the displacements many times, the value of $R$ will trend towards an average, which can be indicated as $\langle R\rangle$. Hence:

$$\langle X^2 \rangle = \langle Y^2 \rangle = \langle Z^2 \rangle = \frac{1}{3}\langle R^2 \tag{2}\rangle$$

There is no bias in any direction, so the values of each of the components are equal. In addition, because $R$ is composed of the sum of three of these components, its average is $1/3$ of each component.
"""

# ╔═╡ 670b6029-7c3d-4f9d-a218-1dbf402667d0
# This is just plotting code, no need to show.
# The plot is similar to the 3D one, but in 2D to make illustration easier.
begin
	if wglmakie
		WGLMakie.activate!()
	else 
		CairoMakie.activate!()
	end
	Random.seed!(0)
	f1 = Figure()
	a1 = Axis(f1[1,1], title="Einstein relation in two dimensions")
	for particle in 1:particles
		x = [0]
		y = [0]
		for step in 1:steps
			# There are four possible directions, choose one
			direction = rand(1:4)
			@match direction begin
				1 => begin
					append!(x, x[end] + 1)
					append!(y, y[end])
				end
				2 => begin
					append!(x, x[end] - 1)
					append!(y, y[end])
				end
				3 => begin
					append!(x, x[end])
					append!(y, y[end] + 1)
				end
				4 => begin
					append!(x, x[end])
					append!(y, y[end] - 1)
				end
				
			end
		end
		lines!(a1, x, y)
		vlines!(a1, [-10, 0, 10], linestyle=:dash)
	end
	f1
end

# ╔═╡ 626913f7-4038-4af0-bf7d-974a85f9cffc
md"""
Adjust the step and particle count to make visuals clearer:

 $step_slider $(steps == 1 ? "step" : "steps")

 $particle_slider $(particles == 1 ? "particle" : "particles").
"""

# ╔═╡ d23e4178-1741-4619-a54c-7d505996ff5e
md"""
Notice that as you increase the number of steps (use more particles to aid in illutration), more particles approach, then cross the planes shown at $x = -10$ and $10$. The number of steps is a function of time, $t$. 

For a plane at $x$, the number of diffusing particles at a time $t + \tau$ is the sum of particles coming from all other planes, where $\tau$ is some time interval. For example, particles that travel from the planes at $X = 10$ and $X = -10$ contribute to the number of particles passing through the plane at $x = 0$. This is true for all values of $X \in \mathbb{Z}$ (All integer values of $X$ from $-\infty$ to $\infty$). 

If the concentration of particles at the plane $x - X$ is $C(x - X, t)$, then the diffusion contribution to the concentration of particles at the $x$-plane ($C(x, t)$) is $C(x - X, t)p(X,\tau)$. $p(X,\tau)$ is the probability that a particle diffuses over a distance $X$ in a time interval $\tau$.

In concrete terms, the plane at $x = -10$ may have $10$ particles passing through it. Over the course of $5$ steps, it will contribute $10\times p(10, 5 \text{ steps})$ particles to the plane at $x = 0$. So, if the probability that a particle travels 10 units of distance over 5 steps is 0.1, the $x = -10$ plane will contribute 1 particle to the $x = 0$ plane.

You can represent the concentration at $x$ as:

$$C(x,t+\tau) = \sum_{X' \in \mathbb{Z}} C(x - X, t)p(X, \tau) \tag{3}$$

where $X'$ is all _other_ planes' contributions to the plane at $x$. The Taylor expansion of concentration yields:

$$C(x,t) + \tau\frac{\delta C}{\delta t} + \cdots = \sum_{X \in \mathbb{Z}} \left[C(x, t) - X\frac{\delta C}{\delta x} + \frac{X^2}{2}\frac{\delta^2C}{\delta x^2} + \cdots\right]p(X, \tau) \tag{4}$$

We can also observe that at any given plane, the probability that a particle will "diffuse" anywhere to another $X$, even if it doesn't move over a given interval ($X = 0$, the total must be equal to unity (or 1). This is true so long as additional particles don't disappear or appear:

$$\sum_X p(X, \tau) = 1 \tag{5}$$


We'll now introduce a useful term, the $n$th moment of $X$:

$$\sum_X X^np(X,\tau) = \langle X^n \rangle \tag{6}$$

Comparing terms in Equation (4) and neglecting terms higher than first order on the left and second order on the right we can conclude:

$$C(x,t) + \overbrace{\tau\frac{\delta C}{\delta t}}^{\text{First order}} + \cdots = \sum_{X \in \mathbb{Z}} \left[C(x, t) \overbrace{- X\frac{\delta C}{\delta x} + X^2\frac{\delta^2C}{\delta x^2}}^{\text{First + second order}} + \cdots\right]p(X, \tau) \tag{7}$$


$$\frac{\delta C}{\delta t} \approx \frac{\langle X\rangle}{\tau} \frac{\delta C}{\delta x} + \frac{\langle X^2\rangle}{2\tau} \frac{\delta^2C}{\delta x^2} \tag{8}$$

The first term on the right $\frac{\langle X\rangle}{\tau} \frac{\delta C}{\delta x}$ can be thought of as a drifting driven by some external force. However, for a true random walk, the average value of $X$ is equally biased to left and right of zero, so $\langle X\rangle = 0$, eliminating this term. This is not the case for $\langle X^2\rangle$. Hence, we can write:

$$\frac{\delta C}{\delta t} = \frac{\langle X^2\rangle}{2\tau} \frac{\delta^2C}{\delta x^2} \tag{9}$$

Note the resemblance to Fick's second law:

$$\frac{\delta C}{\delta t} = D\frac{\delta^2C}{\delta x^2} \tag{10}$$

Looking back to Equation (2), the diffusion constant can be written:

$$D = \frac{\langle X^2\rangle}{2\tau} = \frac{\langle R^2\rangle}{3} \tag{11}$$
"""

# ╔═╡ d1766538-b8d7-458b-951e-01b23a3924cd
md"""
## Random walk

Consider the vector $\vec R$, which is 	composed of the sum of $n$ individual jump vectors $\vec r$. In addition, consider the $x$ component of this vector, $X$, which is composed of the sum of the $x$-component of the jump vectors $\vec r$.

$$\vec R = \sum_{i=1}^n\vec r_i \tag{12}$$


$$X = \sum_{i=1}^n x_i \tag{13}$$

Squaring the magnitudes:

$$|\vec R|^2 = \sum_{i=1}^n|\vec r_i|^2 + 2\sum_{i \ne j} |\vec r_i||\vec r_j| \tag{14}$$


$$X^2 = \sum_{i=1}^n x_i^2 + 2\sum_{i\ne j}^n x_ix_j \tag{15}$$

Looking at Equation (14), for true random jumps, each jump can be treated independently. If we assume each jump takes a distance $d$: $\sum_{i=1}^n |\vec r_i|^2 = nd^2$. In addition, the second term is a sum of a random distribution of values of $|\vec r_i||\vec r_j|$, which can range from $-d^2$ to $d^2$ with an unbiased probability around zero, so this term vanishes. Hence:

$$\langle|\vec R|^2\rangle = \sum_{i=1}^n \langle |r_i|^2\rangle = nd^2 \tag{16}$$

$$\langle X^2\rangle = \sum_{i=1}^n \langle x_i^2\rangle = nd_x^2 \tag{17}$$

We can substitute Equation (16-17) into Equation (11):

$$D = \frac{nd_x^2}{2\tau} = \frac{nd^2}{6\tau}\tag{18}$$

We can further define the number of jumps $n$:

$$n  = \Gamma Z \tau \tag{19}$$

$\Gamma$ is the frequency in which jumps occur, $Z$ is the number of neighboring sites available, and $\tau$ is a time interval. We can define $\Gamma$ using an Arrhenius relation:

$$\Gamma = \Gamma_0\exp\left(-\frac{E_m}{kT}\right) \tag{20}$$

Where $\Gamma_0$ is the number of "attempted" jumps, $E_m$ is the activation energy for a successful jump. $k$ is the Boltzmann constant and $T$ is temperature. The behavior of this function as a function of temperature can be observed in this figure:
"""

# ╔═╡ 51869391-d0e7-408a-9d8d-873afa01976c
# Plotting code
begin
	if wglmakie
		WGLMakie.activate!()
	else 
		CairoMakie.activate!()
	end
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "Temperature (K)", ylabel = "Γ", title="Plot of Γ as a function of temperature")
	ts = range(1, 1000, 100)u"K" # Temperature K
	gamma_0 = 1e14 # jumps/s
	k = BoltzmannConstant
	E_m = 0.1u"eV"
	res = @. E_m / k / ts
	ys = @. gamma_0 * exp(-E_m / k / ts)
	lines!(ax, ts, ys)
	fig
end

# ╔═╡ 3428b15a-6d55-4bcf-b361-00af9d3de667
md"""
The diffusivity can then be expressed as:

$$D = \frac{1}{6}\Gamma Z d^2 \tag{21}$$

Or:

$$D = D_0\exp\left(-\frac{E_m}{k_T}\right) \text{ with } D_0 = \frac{1}{6}\Gamma_0 Z d^2 \tag{22}$$
"""

# ╔═╡ 8518ab0e-173f-4ee2-a0a8-7e3eeaf4fcc6
md"""
## Solid state diffusion correlation

In solids, jumps may be correlated. For example, a jump leaves a vanacy, and the particle may be biased to jumping back to the vacancy it had left behind. We can represent this correlation as follows:

$$\langle|\vec R_c|\rangle^2 = nd^2 + d^2 2 \sum_{i,j}\cos \theta_{i,j} \tag{23}$$

 $\theta_{i,j}$ is the angle between the jump vectors $\langle r_i \rangle$ and $\langle r_j \rangle$. The subscript $c$ refers to correlated diffusion. This is the case in which $\sum_{i,j}\cos\theta_{i,j} \ne 0$.

We can define a correlation factor $f$ that relates the ratio of $\langle |\vec R_c|^2\rangle$ to $\langle |\vec R|^2\rangle$:

$$f = \frac{\langle |\vec R_c|^2\rangle}{\langle |\vec R|^2\rangle} = 1 + \frac{2\sum_{i,j}\cos \theta_{i,j}}{n} \tag{24}$$. A more general diffusivity can then be written:

$$D = \frac{1}{6}\Gamma Z f d^2 \tag{25}$$

### Interstitial diffusion

Interstitials (atoms stuck between lattice sites in a crystal) diffuse randomly. They do not have a correlation between jumps. When an interstitial jumps to a new spot, it is not restricted from moving to any neighboring spot after the jump. For interstitial diffusion, $f = 1$.


"""

# ╔═╡ ba5df24e-436b-4269-b355-8d52e7a00be2
# Plotting code
begin
	function interstitial_diffusion()
		CairoMakie.activate!()
		fig = Figure()
		axis = Axis(fig[1,1], title="Interstitial Diffusion")
		xlims!(axis, -0.25, 4.25)
		ylims!(axis, -0.25, 4.25)
		# Generate grid coordinates
		x1 = 0:4
		y = 0:4
		
		# Create all combinations using vec and repeat
		x = vec([i for i in x1, j in y])
		y = vec([j for i in x1, j in y])
	
		record(fig, "interstitial_diffusion.gif", 1:100, framerate=30) do i
			empty!(axis)
			scatter!(axis,  x, y, markersize=50)
			
			offset = i > 50 ? ((100 - i) / 50) : i / 50
			scatter!(axis, [1.5 + offset], [1.5], markersize=50)
		end
		LocalResource("interstitial_diffusion.gif")
	end
	interstitial_diffusion()
end

# ╔═╡ 63cd3890-db15-4239-9fae-0e8565f8a7ac
md"""
### Vacancy diffusion

Vacancies can only diffuse by swapping places with an adjacent atom. If an atom jumps to occupy a vacancy, there is a high likelihood that the atom will return to the original site. In the case that the atom returns to the original site, $\cos \theta_{i,i+1} = -1$. 

Imagine in a 3D lattice, a vacancy diffuses by swapping with a "trace" atom. Consider the yellow atom in the animation below as a "trace" atom. There is $Z$ neighboring sites available for the vacancy to diffuse to. There is a $1/Z$ chance that the vacancy swaps again with the trace atom. Otherwise, there is a $(Z-1)/Z$ chance that the vacancy will diffuse elsewhere, and the trace atom remains displaced.

The expected value of $\langle \cos \theta \rangle$ for adjacent jumps is:

$$\langle \cos \theta \rangle = \frac{1}{Z}(-1) + \frac{(Z-1)}{Z}(0) = -\frac{1}{Z} \tag{26}$$

And for large $n$ and only considering adjacent jumps, $\langle \cos\theta\rangle$ can also be approximated:

$$\langle \cos \theta \rangle \approx \frac{\sum_{i,j}^n\cos\theta_{i,j}}{n} \tag{27}$$

Plugging into Equation (24), the vacancy diffusion mechanism can be represented:

$$f \approx 1 - 2\frac{1}{Z} \tag{28}$$

Note that this simplified model does not account for vacancies returning to their original place after multiple steps.
"""

# ╔═╡ 891ace0a-83b4-4f85-b571-ad19969aea22
# Plotting code
begin
	function vacancy_diffusion()
	CairoMakie.activate!()
	fig = Figure()
	axis = Axis(fig[1,1], title="Vacancy Diffusion (1/Z) case")
	
	xlims!(axis, -0.25, 4.25)
	ylims!(axis, -0.25, 4.25)

	x1 = 0:4
	y = 0:4
	
	# Create all combinations using vec and repeat
	x = vec([i for i in x1, j in y])
	y = vec([j for i in x1, j in y])

	x = deleteat!(x, [13, 14])
	y = deleteat!(y, [13, 14])
	
	record(fig, "vacancy_diffusion.gif", 1:100, framerate=30) do i
		empty!(axis)
		scatter!(axis,  x, y, markersize=50)
		
		offset = i > 50 ? ((100 - i) / 50) : i / 50
		x2 = [2 + offset]
		scatter!(axis, x2, [2], markersize=50)
		text!(axis, x2, [2.2], text = "Trace Atom")
		text!(axis, [3 - offset], [1.8], text = "Vacancy")
		
	end
	LocalResource("vacancy_diffusion.gif")
	end
	vacancy_diffusion()
end

# ╔═╡ c1a1aa76-1b8d-4abf-a697-83d9ffab3875
# Plotting code
begin
	function vacancy_diffusion_alt()
		CairoMakie.activate!()
		fig = Figure()
		axis = Axis(fig[1,1], title="Vacancy Diffusion (Z-1)/Z case")
	
		xlims!(axis, -0.25, 4.25)
		ylims!(axis, -0.25, 4.25)
			
		x1 = 0:4
		y = 0:4
			
		# Create all combinations using vec and repeat
		x = vec([i for i in x1, j in y])
		y = vec([j for i in x1, j in y])
	
		x = deleteat!(x, [13, 14, 18])
		y = deleteat!(y, [13, 14, 18])
	
		record(fig, "vacancy_diffusion_alt.gif", 1:150, framerate=30) do i
			empty!(axis)
			scatter!(axis,  x, y, markersize=50)
	
			if i < 50
				offset = i / 50
				x2 = [2 + offset]
				scatter!(axis, x2, [2], markersize=50)
				text!(axis, x2, [2.2], text = "Trace Atom")
				text!(axis, [3 - offset], [1.8], text = "Vacancy")
				scatter!(axis, [2], [3], markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
			elseif i < 100
				offset = (i - 50) / 50
				y2 = [3 - offset]
				scatter!(axis, [3], [2], markersize=50)
				text!(axis, [3], [2.2], text = "Trace Atom")
				text!(axis, [2], [1.8 + offset], text = "Vacancy")
				scatter!(axis, [2], y2, markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
			else
				scatter!(axis, [3], [2], markersize=50)
				text!(axis, [3], [2.2], text = "Trace Atom")
				text!(axis, [2], [2.8], text = "Vacancy")
				scatter!(axis, [2], [2], markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
			end
		end
		LocalResource("vacancy_diffusion_alt.gif")
		end
	vacancy_diffusion_alt()
end


# ╔═╡ 8f6d9c95-0cfa-41c1-ae3c-c5e7ec6d8fd4
md"""
## Macroscopic diffusivity

The previously described diffusivity, $D_m = \frac{1}{6}\Gamma Z fd^2$ is the _microscopic diffusivity_. This is not the same as the _macroscopic diffusivity_, which is measured experimentally. The macroscopic diffusivity accounts for a mobile atom diffusing, but it also accounts for stationary atoms becoming diffusing atoms, and diffusing atoms becoming stationary.

We can write the macroscopic diffusivity:

$$D = \lambda^2g = g\left[\sqrt{\frac{D_m}{r}}\right]^2 = \frac{1}{6}\Gamma Z f d^2 \frac{g}{r} \tag{29}$$

where $\lambda^2 = \frac{D_m}{r}$ is the mean square distance an atom diffuses while it is mobile, $r$ is the recombination rate, or the rate at which mobile atoms becomes stationary, and $g$ is the generation rate, or the rate at which stationary atoms become mobile.
"""

# ╔═╡ afc0732f-f49a-4b27-ade6-88e0936ccd5a
# Plotting code
begin
	function mobile_stationary_plot()
		CairoMakie.activate!()
		fig = Figure()
		axis = Axis(fig[1,1], title="Mobile/stationary atoms")

		xlims!(axis, -0.25, 4.25)
		ylims!(axis, -0.25, 4.25)

		x1 = 0:4
		y = 0:4
		
		# Create all combinations using vec and repeat
		x = vec([i for i in x1, j in y])
		y = vec([j for i in x1, j in y])
	
		deleteat!(x, 13)
		deleteat!(y, 13)
	
		record(fig, "mobile_stationary_diffusion.gif", 1:150, framerate=30) do i
			empty!(axis)
	
			if i < 50
				offset =  i / 100
	
				x2 = [2 + offset]
				y2 = [2 + offset]
				
				scatter!(axis, x, y, markersize=50)
				scatter!(axis, x2, y2, markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
				text!(axis, [2 + offset], [2.2 + offset], text = "Stationary Atom")
				scatter!(axis, [1.5], [1.5], markersize=50)
				text!(axis, [1.5], [1.7], text = "Mobile Atom")
			elseif i < 100
				offset =  (i - 50) / 100
				
				y2 = [1.5 + offset]
				x2 = [1.5 + offset]
				scatter!(axis, x, y, markersize=50)
				scatter!(axis, [2.5], [2.5], markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
				scatter!(axis,  x2, y2, markersize=50)
				text!(axis, [2.5], [2.7], text = "Stationary Atom")
				text!(axis, [1.5 + offset], [1.7 + offset], text = "Mobile Atom")
			else
				scatter!(axis, x, y, markersize=50)
				scatter!(axis, [2.5], [2.5], markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
				scatter!(axis,  [2], [2], markersize=50)
				text!(axis, [2.5], [2.7], text = "Mobile Atom")
				text!(axis, [2], [2.2], text = "Stationary Atom")
			end
		end
		LocalResource("mobile_stationary_diffusion.gif")
		end
	mobile_stationary_plot()
end

# ╔═╡ ab388d5e-651d-4962-9b7e-1dfdf39b35fe
md"""
### Diffusivity of point defects

Diffusing point defects themselves can be treated as always being mobile (neglecting conversion to lattice atoms), and the activation and deactivation rates can be taken to be unity:

$$D_i = \frac{1}{6}\Gamma Z f d^2 \tag{30}$$

## Diffusivity of defect assisted atomic diffusion

A trace atom may diffuse via vacancy diffusion, in which case the generation rate $g$ is determined by the availability of a neighboring empty site (vacancy). The likelyhood that any given site is vacant can be derived using physical chemistry:

$$p_V = g_V\exp\left(-\frac{G_V^F}{kT}\right)=  g_v\exp\left(\frac{S_V^F}{k}\right)\exp\left(-\frac{H_V^F}{kT}\right) \tag{31}$$

where $G_V^F$ is the Gibb's free energy of vacancy formation, $S_V^F$ is the entropy of vacancy formation, $H_V^F$ is the enthalpy of vacancy formation, $k$ is the Boltzmann constant, $T$ temperature, and $g_V$ is the geometric factor. This accounts for multiple states For vacancy formation in a monoatomic crystal, $g_v = 1$. 

The derivation of this expression will not be given here, but can be found online:

1. [engineeringenotes.com](https://www.engineeringenotes.com/metallurgy/defects-metallurgy/equilibrium-concentration-of-vacancies-in-crystals-metallurgy/43153)
2. [morfli.com](https://www.morfli.com/lesson-page/1543)
3. [dtrinkle.matse.illinois.edu](https://dtrinkle.matse.illinois.edu/MatSE584/kap_2/exercise/s2_1_4.html)

The generation rate is the product of the number of available sites $Z$ and $p_V$, giving:

$$g = Zp_v\tag{32}$$
"""

# ╔═╡ 7c249d66-6ab2-4161-977e-abc9f719073a
# Plotting code
begin
	function vacancy_formation()
		CairoMakie.activate!()
		fig = Figure()
		axis = Axis(fig[1,1], title="Vacancy formation")

		xlims!(axis, -0.25, 4.25)
		ylims!(axis, -0.25, 4.25)
			
		x1 = 0:4
		y = 0:4
		
		# Create all combinations using vec and repeat
		x = vec([i for i in x1, j in y])
		y = vec([j for i in x1, j in y])
	
		deleteat!(x, 13)
		deleteat!(y, 13)
	
		record(fig, "vacancy_formation.gif", 1:100, framerate=30) do i
			empty!(axis)
	
			if i < 50
				offset =  i / 100
	
				x2 = [2 + offset]
				y2 = [2 + offset]
				
				scatter!(axis, x, y, markersize=50)
				scatter!(axis, x2, y2, markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
				text!(axis, [2 + offset], [2.2 + offset], text = "Thermally excited atom")
			else
				scatter!(axis, x, y, markersize=50)
				scatter!(axis, [2.5], [2.5], markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
				text!(axis, [2.3], [2.7], text = "Interstitial")
				text!(axis, [1], [2.2], text = "Potentially Mobile")
				text!(axis, [2], [1.8], text = "Vacancy")
			end
		end
		LocalResource("vacancy_formation.gif")
	end
	vacancy_formation()
end

# ╔═╡ b0c71984-34cc-4496-a1a9-e2a01ac1dea7
md"""
When a vacancy diffuses away, an atom that was mobile may cease to be mobile. The recombation rate is given by the rate at which a vacancy will migrate:

$$r = r_0\exp\left(-\frac{G_V^M}{kT}\right) = r_0\exp\left(\frac{S_V^M}{k}\right)\exp\left(-\frac{H_V^M}{kT}\right) \tag{33}$$

where $G_V^M$ is the Gibb's free energy of vacancy formation, $S_V^M$ is the entropy of vacancy formation, $H_V^M$ is the enthalpy of vacancy formation, and $r_0$ is the rate at which a migration is attempted.
"""

# ╔═╡ 0188a13e-0b79-496c-992b-74227047529f
# Plotting code
begin
	function vacancy_migration()
		CairoMakie.activate!()
		fig = Figure()
		axis = Axis(fig[1,1], title="Vacancy migration")

		xlims!(axis, -0.25, 4.25)
		ylims!(axis, -0.25, 4.25)
	
		x1 = 0:4
		y = 0:4
		
		# Create all combinations using vec and repeat
		x = vec([i for i in x1, j in y])
		y = vec([j for i in x1, j in y])
	
		deleteat!(x, [13,14])
		deleteat!(y, [13,14])
	
		record(fig, "vacancy_migration.gif", 1:100, framerate=30) do i
			empty!(axis)
	
			if i < 50
				offset =  i / 50
	
				x_alt = [3 - offset]
				y_alt = [2]
				
				scatter!(axis, x, y, markersize=50)
				scatter!(axis, x_alt, y_alt, markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
				text!(axis, [2 + offset], [2.2], text = "Vacancy")
				text!(axis, [1], [2.2], text = "Mobile")
			else
				scatter!(axis, x, y, markersize=50)
				scatter!(axis, [2], [2], markersize=50, color = 1, colormap = :tab10, colorrange = (1, 10))
				text!(axis, [3], [2.2], text = "Vacancy")
				text!(axis, [1], [2.2], text = "Not Mobile")
			end
		end
		LocalResource("vacancy_migration.gif")
		end
	vacancy_migration()
end

# ╔═╡ dcd81bf5-c7e4-4c93-a7f4-6657d638e0b4
md"""
## Nernst-Einstein equation

In thermodynamic equilibrium, the distribution of impurities and defects follows the Boltzmann distribution:

$$C(x) = C_0\exp\left(-\frac{E(x)}{kT}\right) \tag{34}$$

where $C(x)$ is the concentration of impurities or defects as a function of position, $C_0$ is the pre-exponential factor, $E(x)$ is energy, $k$ is the Boltzmann constant, and $T$ is the temperature.

If energy is a function of position, the Boltzmann distribution will create a higher concentration of particles in the low-energy regions. As a result, diffusion will create a flux from high concentration to low concentration.

At equilibrium, there must be no net flux. As particles diffuse away from the high concentration areas, and thus energy is carried away, a drift velocity in the opposite direction must be present to balance the particles lost. A _force_ pushes particles back into the high concentration areas to balance the flux. Fick's law can then be written:

$$J = -D\frac{\delta C}{\delta x} + \bar v C = 0 \tag{35}$$
$$\bar v = \frac{D}{C}\frac{\delta C}{\delta x} \tag{36}$$

The drift velocity $\bar v$ is established by $E(x)$. Substituting Equation (34) into Equation (36) gives:

$$\bar v = \frac{D\left(-\frac{\delta E}{\delta x}\right)}{kT} = \frac{DF}{kT} \tag{37}$$

Where $F$ is the force acting on a particle due to the energy/potential gradient. This is the Nernst-Einstein equation.
"""

# ╔═╡ 34c17de5-576a-40ba-b5a4-afa9ee2d8aef
md"""
## Defect Trapping

There are a few ways in which defects can be trapped. When two defects meet, they may recombine, or agglomerate into a larger defect cluster. The rate at which trapping occurs depends on the configuration of a defect cluster. The rate of trapping will be different for a spherical vs. planar vs. rod-like defect cluster, for example.

When modeling trapping, two assumptions are made:

1. Trapping occurs when the separation distance between two defects are less than a critical radius.
2. Within the critical radius, the defect concentration decreases linearly and reaches zero when the radius $r=0$.

In the case of a spherical defect cluster $A$ with a trapping radius $r$, the trapping flux at $r$ is determined by defect diffusivity and concentration gradient at $r$ of a defect $B$. Assuming linear depletion, the trapping flux is given:

$$J_{B\rightarrow A} = -D_B\frac{\delta C_B}{\delta x} = -D_B\frac{C_B(r) - 0}{r} \tag{38}$$

To find the total trapping rate of $B$, one has to multiply by the surface area of available trapping sites (with concentration $C_A$) on the spherical defect cluster by the flux:

$$\left[\frac{\Delta C_B}{\Delta t}\right]_{\text{spherical}} = C_A4\pi r^2 J_{B\rightarrow A} = (4\pi rD_BC_B)C_A \tag{39}$$

$4\pi r$ represents the effective trapping area of a single defect $A$ with a radius $r$.

"""

# ╔═╡ 5b23e294-6bf6-439c-84ff-34cf99d31e53
# Plotting code
begin
function defect_tracking()
	CairoMakie.activate!()
	fig = Figure()
	axis = Axis(fig[1,1], title="Defect trapping")
	
	xlims!(axis, -0.25, 4.25)
	ylims!(axis, -0.25, 4.25)

	x1 = 1.5:0.25:2.5
	y = 0:0.5:4
	
	# Create all combinations using vec and repeat
	x = vec([i for i in x1, j in y])
	y = vec([j for i in x1, j in y])

	deleteat!(x, 25);
	deleteat!(y, 25);

	record(fig, "defect_trapping.gif", 1:100, framerate=30) do i
		empty!(axis)

		scatter!(axis, x, y, markersize=50)
		text!(axis, [2.75], [4], text = "Rod-Like Defect Cluster")
		
		if i < 50
			offset =  i / 50
			 
			x_alt = [3.5 - offset]
			y_alt = [2]

			scatter!(axis, x_alt, y_alt, markersize=50)
			lines!(axis, [2.25, x_alt[1]], [2, 2], color=:black)
			text!(axis, [3.5 - offset], [2.2], text = "Defect")
			text!(axis, [3.25 - offset], [1.8], text = "Separation = $(round(sigdigits=3, 1.5 - offset))")
		else
			scatter!(axis, [2.5], [2], markersize=50)
			text!(axis, [2.5], [2.2], text = "Defect", color=:black)
			lines!(axis, [2.25, 2.5], [2, 2], color=:black)
			text!(axis, [2.25], [1.8], text = "Separation = 0.5")
		end
	end
	LocalResource("defect_trapping.gif")
end
	defect_tracking()
end

# ╔═╡ c1feebf7-4384-4196-bede-8f0662f301c4
md"""
For a 2D dislocation loop (to be treated as a circle), the trapping flux $J_{B\rightarrow A}$ takes the same form as in Equation (38). But the effective trapping area consists of the surface area of a circle instead with both sides available to absorb defects, which is $2\times \pi R^2$, where $R$ is the loop radius:

$$\left[\frac{\Delta C_B}{\Delta t}\right]_{\text{loop}} = C_A(2\pi R^2)D_B\frac{C_B}{r} = \left(2\pi \frac{R^2}{r}D_BC_B\right)C_A \tag{40}$$

For a 1D rod-like defect cluster, the ends are ignored and the surface area is $(\pi RL$, where $R$ is the rod radius, and $L$ is the length of the rod.

$$\left[\frac{\Delta C_B}{\Delta t}\right]_{\text{rod}} = C_A(2\pi RL)D_B\frac{C_B}{r} = \left(2\pi\frac{RL}{r}D_BC_B\right) \tag{41}$$

"""

# ╔═╡ 9ae89154-24de-42b6-81ae-ebb07d6a87ea
# Plots different types of defect clusters
begin
	function defect_cluster_types()
		if wglmakie
			WGLMakie.activate!()
		else 
			CairoMakie.activate!()
		end
		fig = Figure()
		ax1 = Axis3(fig[1, 1], title="Spherical Defect")
		ax2 = Axis3(fig[1, 2], title="Loop (circular) Defect")
		ax3 = Axis3(fig[1, 3], title = "Rod-like defect")
		xlims!(ax3, -1, 1)
		ylims!(ax3, -1, 1)
		ax4 = Axis(fig[2, 1:3], xlabel = "Radius", ylabel = "Concentration", title = "Defect concentration vs. radius")

		sphere = Sphere(Point3f(0, 0, 0), 1.0f0)
		mesh!(ax1, sphere, color = :dodgerblue)

		circle = Circle(Point2f(0,0), 1.0f0)
		mesh!(ax2, circle, color = :dodgerblue)

		rod = Cylinder(Point3f(0, 0, -2), Point3f(0, 0, 2), 0.3f0)
		mesh!(ax3, rod, color = :dodgerblue)

		x = [1,2,5]
		y = [0,1,1]
		lines!(ax4, x, y)
		vlines!(ax4, 2, linestyle=:dash)
		
		fig
	end
	defect_cluster_types()
end

# ╔═╡ 9ca2fdcf-5909-4c52-b470-080add85212a
md"""
In each of these three cases, linear depletion is assumed. This approximation does not consider the conservation of particles. In reality, for trapping by spherical and rod-like defect clusters, the flux is meant to increase as radius decreases, as the available surface area decreases. But this effect is negligible, and the estimation performs well enough without accounting for it.

In the above derivation, it was assumed that the defect cluster was immobile. In the case of two mobile point defects (interstitials combining with vacanceis), the expression for trapping rate becomes:

$$\left[\frac{\Delta C_I}{\Delta t}\right]_{I-V} = \left[\frac{\Delta C_I}{\Delta t}\right]_{I-V} = 4\pi r^2(D_I + D_V)C_IC_V \tag{41}$$

Where $I$ and $V$ are vacancies and interstitials, respetively. You are encouraged to attempt to make this derivation yourself.

The trapping radius is generally taken to be the value of the lattice spacing. In addition, it is assumed that there is no energy barrier, and that a defect will freely become trapped once inside the trapping radius. In the case that a trapping reaction does require energy, the expression of the trapping radius can be augmented:

$$r = r_0\exp\left(-\frac{E_T}{kT}\right)\tag{42}$$

where $r_0$ is the energy-free trapping radius, $E_T$ is the trapping activation energy, $k$ is the Boltzmann constant, and $T$ is temperature.
"""
# I don't really know what this means.

# ╔═╡ dc173930-d763-4e53-baca-be8779f41614
begin
	function lattice_spacing()
		if wglmakie
			WGLMakie.activate!()
		else 
			CairoMakie.activate!()
		end

		fig = Figure()
		ax = Axis3(fig[1,1], title = "Lattice spacing (BCC)")

		x = [0, 0, 0, 0, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 1]
		y = [0, 0, 1, 1, 0, 0, 1, 1, 0.5, 0.5, 0.5, 0, 1, 0.5, 0.5]
		z = [0, 1, 0, 1, 0, 1, 0, 1, 0.5, 0, 1, 0.5, 0.5, 0.5, 0.5]
		scatter!(ax, x, y, z, markersize=50)
		lines!(ax, [0,0,0], [0,1,0], linestyle=:dash)
		text!(ax, 0, 0.5, 0, text="r = Lattice spacing")
		
		
		fig
	end
	lattice_spacing()
end

# ╔═╡ 591e64a6-56d6-4014-8db6-d90274f5216f
md"""
## Defect interactions

Excluding large defect clusters, the concentration of interstitials and vacancies can be modeled:

$$\frac{\delta C_I}{\delta t} = \nabla\cdot D_I\nabla C_I + K_0 - K_{I-V}C_IC_V - K_{I-T}C_IC_T  \tag{43}$$
$$\frac{\delta C_V}{\delta t} = \nabla\cdot D_V\nabla C_V + K_0 - K_{I-V}C_IC_V - K_{V-T}C_VC_T  \tag{44}$$

Where $K_0$ is the generation rate of point defects, $K_{I-V}[=4\pi r(D_I + D_V)]$ is the defect recombination rate, $D_{I,V}$ and $C_{I,V}$ are diffusivity and concentration of interstitials $(I)$ and vacancies $(V)$.
"""

# ╔═╡ 2322c0d0-3fa4-47fd-8c3c-905af02461d5
begin
	function c_over_time_graph()
		CairoMakie.activate!()
		

		fig = Figure()
		ax = Axis(fig[1,1], xlabel = L"\log t", ylabel = L"\log C", title = "Vacancy/Interstitial concentration over time")
		xlims!(ax, 0, 4.25)
		ylims!(ax, 0, 2)
	
		lines!(ax, [0,1], [0,1], color=:black)
		vlines!(ax, 1, ymin = 0, ymax = 0.5, linestyle = :dash, color = :black)
		text!(ax, 0.2, 1, text="low temperature\nlow sink density\nτ₁ < τ₂")
		text!(ax, 0.2, 0.85, text=L"C_I = C_v = K_0t")
		text!(ax, 0.95, 0.3, text=L"(K_0K_{I-V})^{-1/2}", rotation=π/2)

		lines!(ax, [1,2], [1,1], color=:black)
		text!(ax, 1, 1, text=L"C_I = C_V = (K_0/K_{I-V})^{1/2}")
		vlines!(ax, 2, ymin = 0, ymax = 0.5, linestyle = :dash, color = :black)
		text!(ax, 1.95, 0.3, text=L"(K_{I_T}C_T)^{-1}", rotation=π/2)
		
		lines!(ax, [2,3], [1,1.5], color=:black)
		lines!(ax, [2,3], [1,0.5], color=:black)
		vlines!(ax, 3, ymin = 0, ymax = 0.75, linestyle = :dash, color = :black)
		text!(ax, 1.8, 1.5, text=L"C_V = (K_0K_{I-S}C_St / K_{I-V})^{1/2}")
		text!(ax, 2.0, 0.3, text=L"C_I = \left(\frac{K_0}{K_{I-V}L_{I-T}C_Tt}\right)^{1/2}")
		text!(ax, 3, 0.75, text=L"(K_{V-T}C_T)^{-1}", rotation = π/2)

		lines!(ax, [3,4], [1.5,1.5], color=:black)
		lines!(ax, [3,4], [0.5,0.5], color=:black)
		text!(ax, 3.1, 1.5, text=L"C_V = \left(\frac{K_0K_{I-T}}{K_{I-V}K_{V-T}}\right)^{1/2}")
		text!(ax, 3.1, 0.5, text=L"C_I = \left(\frac{K_0K_{V-T}}{K_{I-V}K_{I-T}}\right)^{1/2}")
		text!(ax, 0, -0.5, text=L"C_I = \left(\frac{K_0K_{V-T}}{K_{I-V}K_{I-T}}\right)^{1/2}")

		text!(ax, 0.1, 1.7, text= "buildup without\nreaction")
		text!(ax, 1.1, 1.7, text= "mutual\nrecombination\ndominates")
		text!(ax, 1.1, 1.7, text= "mutual\nrecombination\ndominates")
		text!(ax, 2.1, 1.7, text= "sinks contribute\nto interstitial\nannihilation")
		text!(ax, 3, 1.8, text= "sinks also contribute\n to vacancy annihilation")
		
		fig
	end
	c_over_time_graph()
end

# ╔═╡ 16d9d932-ddf4-4486-8324-ae6dcdf86e18
md"""
## Rate theory for defect clustering

The above  approximations are adequate for predicting damage without defect clustering, which is valid for low temperature or low defect concentration. When irradition takes place over a long period of time, defects begin to agglomerate into larger defect clusters. This is especially prevelant with longer annealing time and high temperature annealing. 

These clusters can be accounted for with the previous equations by adding terms which account for the absorption and release of defects.

$$\frac{\delta C_I}{\delta t} = \nabla\cdot D_I\nabla C_I + K_0 - K_{I-V}C_IC_V - K_{I-T}C_IC_T - \underbrace{\sum_{n\ge 2}K_{I-I_n}C_IC_{I_n} + \sum_{n\ge 2}R_{I_n}C_{I_n}\tag{45}}_{\text{Interstitial cluster absorption/emission}}$$
$$\frac{\delta C_V}{\delta t} = \nabla\cdot D_V\nabla C_V + K_0 - K_{I-V}C_IC_V - K_{V-T}C_VC_T - \underbrace{\sum_{n\ge 2}K_{V-V_n}C_VC_{V_n} + \sum_{n\ge 2}R_{V_n}C_{V_n}}_{\text{Vacancy cluster absorption/emission}} \tag{46}$$

where $K_{I-I_n}, K_{V-V_n}$ is the rate at which interstitials/vacancies are absorbed into the cluster of size $n$, $C_{I_n}, C_{V_n}$ is the concentration of interstitials/vacancies in the cluster, $R_{I_n}, R_{V_n}$ is the rate at which defect clusters decay (by emitting interstitials/vacancies)

Defects will attempt to escape defect clusters, but must have a threshold energy high enough to disaccociate and migrate away from the defect cluster. The required energy is a combiatnion of the diassociation energy and the migration energy. The Boltzmann distribution determines the fraction of defects that can escape a cluster.
"""

# ╔═╡ 4365438f-afa6-4478-9d86-6c82545dc84e
begin
	function energy_graph()
		CairoMakie.activate!()
		
		fig = Figure()
		ax = Axis(fig[1,1], xlabel = "Configuration", ylabel = "Energy", title = "Defect and dissociation energy ")
		hidedecorations!(ax, label=false,)
		
		# xlims!(ax, 0, 4.25)
		ylims!(ax, 0, 1.5)

		lines!(ax, [0,1], [0,1], color=:black)
		x = 1:0.01:3
		y = @. -(sin((x + 1) * 11))/8 .+ 1
		lines!(ax, x, y, color=:black)
		
		x = 3:0.01:4
		y = @. ((x-3.5)*10)^2/50 + 0.5
		lines!(ax, x, y, color=:black)
		
		x = 4:0.01:6
		y = @. (cos((x + 1) * 11))/8 .+ 1
		lines!(ax, x, y, color=:black)

		x = 6:0.01:7
		y = @. ((x-6.5)*6)^2/10+0.11
		lines!(ax, x, y, color=:black)

		x = 7:0.01:8
		y = @. (sin((x + 1) * 11))/8 .+ 1
		lines!(ax, x, y, color=:black)

		annotation!(ax, 180, 0, 0.3, 1,
		    text = " ",
		    path = Ann.Paths.Arc(height = 1),
		    style = Ann.Styles.LineArrow()
		)
		annotation!(ax, -180, 1, 6.4, 1,
		    text = " ",
		    path = Ann.Paths.Arc(height = 1),
		    style = Ann.Styles.LineArrow()
		)
		annotation!(ax, 0, 110, 2, 0.47,
		    text = " ",
		    path = Ann.Paths.Line(),
		    style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(), tail = Ann.Arrows.Head())
		)
		
		annotation!(ax, 200, 0, 1.2,  0.5,
		    text = " ",
		    path = Ann.Paths.Line()
		)
		annotation!(ax, 100, 0, 1.2,  0.87,
		    text = " ",
		    path = Ann.Paths.Line()
		)
		text!(ax, 1.6, 0.65, text = L"E_b")
		text!(ax, 1.3, 0.35, text = "Binding energy")
		
		annotation!(ax, -40, -40, 3, 0.25,
		    text = "Formation\nenergy of defect",
		    path = Ann.Paths.Line(),
		    style = Ann.Styles.LineArrow(head = Ann.Arrows.Head())
		)

		annotation!(ax, 0, 138, 3.5, -0.03,
		    text = " ",
		    path = Ann.Paths.Line(),
		    style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(), tail = Ann.Arrows.Head()))
		text!(ax, 3.0, 0.25, text = L"E_{fc}")
		annotation!(ax, 0, 165, 4.14, 0.47,
		    text = " ",
		    path = Ann.Paths.Line(),
		    style = Ann.Styles.LineArrow(head = Ann.Arrows.Head())
		)
		text!(ax, 4.3, 0.6, text = L"E_{\text{diss}}")

		annotation!(ax, 20, -50, 4.5, 0.62,
		    text = "Dissociation energy\nfor I release",
		    path = Ann.Paths.Line(),
		    style = Ann.Styles.LineArrow(head = Ann.Arrows.Head())
		)
		
		fig
	end
	energy_graph()
end

# ╔═╡ f5ef096e-de1f-4299-a8d0-4f9d757f96a2
md"""
We can write equations for each defect cluster of size $n$ as well:

$$\frac{\delta C_{I_n}}{\delta t} = \nabla\cdot D_{I_n}\nabla C_{I_n} + K_{I-I_{n-1}}C_IC_{I_{n-1}} - K_{I-I_n}C_IC_{I_n}	 - R_{I_n}C_{I_n} +  R_{I_{n+1}}C_{I_{n+1}} \tag{47}$$
$$\frac{\delta C_{V_n}}{\delta t} = \nabla\cdot D_{V_n}\nabla C_{V_n} + K_{V-V_{n-1}}C_VC_{V_{n-1}} - K_{V-V_n}C_VC_{V_n} - R_{V_n}C_{V_n} +  R_{V_{n+1}}C_{V_{n+1}} \tag{48}$$

Here, the first term is the diffusion term, where $D_{I_n}, D_{V_n}$ is the diffusivity of a $n$ defect cluster, and $C_{I_n}, C_{V_n}$ is the concentration of a cluster with $n$ defects. $K_{I-I_{n-1}}, K_{V-V_{n-1}}$ is the rate at which defects are trapped by clusters size $n-1$, causing them to become clusters of size $n$. $K_{I-I_{n}}, K_{V-V_{n}}$ is the rate at which clusters of size $n$ capture defects and become size $n + 1$. $R_{I_n}, R_{V_n}$ is the rate at which defects are emmited (decay) from defect clusters of size $n$, and  $R_{I_{n+1}}, R_{V_{n+1}}$ is the rate at which defects are emitted from defect clusters of size $n + 1$.
 
The diffusion term (first term on the right) is generally negligible for defect clusters. They are largely immobile ($D_{I_n} = D_{V_n} = 0$), with the exception of di-interstitials, and di-vacancies, which are thought to be mobile in some metals.

The trapping and decay terms can be written:

$$K_{I - I_n} = 4\pi r_n(D_I + D_{I_n}) \tag{49}$$
$$R_{I_n} = \Gamma_0Zn\exp\left(-\frac{E_{\text{diss}}(n)}{kT}\right) \tag{50}$$

with analogous equations for vacancies.

Where $\Gamma_0$ is the number of jump attempts along a specific direction,  $Z$ is the number of available jump sites, and $E_{\text{diss}}(n)$ is the dissassociation energy for a defect cluster of size $n$.

"""

# ╔═╡ a1b55842-ef12-44f0-a686-4976cd55765a
# Plotting code
begin
function defect_trapping_emission()
	CairoMakie.activate!()
	fig = Figure()
	ax1 = Axis(fig[1,1], title="Defect trapping/emission")

	xlims!.(ax1, 0, 10)
	ylims!.(ax1, 0, 10)
	
	x = [4,6,4,5,6,4,6]
	y = [4,4,5,5,5,6,6]

	record(fig, "defect_trapping_emission.gif", 1:100, framerate=30) do i
		empty!(ax1)
		scatter!(ax1, x, y, markersize=50)
		
		if i < 50
			offset = i / 50
			scatter!(ax1, [5], [6 + offset], markersize=50)
			annotation!(ax1, 5, 7 + offset, text = L"\text{Escaping Defect} R_{I_n}")
			scatter!(ax1, [5], [3 + offset], markersize=50)
			annotation!(ax1, 5, 2 + offset, text = L"\text{Trapped Defect} K_{I-I_n}")
		else
			scatter!(ax1, [5], [7], markersize=50)
			annotation!(ax1, 5, 8, text = L"\text{Escaping Defect} R_{I_n}")
			scatter!(ax1, [5], [4], markersize=50)
			annotation!(ax1, 5, 3, text = L"\text{Trapped Defect} K_{I-I_n}")
		end
		
		
	end
	LocalResource("defect_trapping_emission.gif")
end
	defect_trapping_emission()
end

# ╔═╡ 7b0f6bf6-eb66-4f1a-bc5c-71cfbe2b0c90
md"""
The interstitial diffusivity can be written:

$$D_I = \frac{1}{6}\Gamma_0Z d^2 \exp\left(-\frac{E_m}{kT}\right)  = D_0\exp\left(-\frac{E_m}{kT}\right) \tag{51}$$

$$E_{\text{diss}}(n) = E_{\text{bind}}(n) + E_m\tag{52}$$

where $E_\text{bind}(n)$ is the binding energy of a defect to the defect cluster of size $n$, and $E_m$ is the migration energy. Recall Equation (22) with $D_0 = \frac{1}{6}\Gamma_0Zd^2$.

The defect emission rate can be written:

$$R_{I_n} = 6\frac{D_I}{d^2}n\exp\left[-\frac{E_{\text{bind}}(n)}{kT}\right] \tag{53}$$

Or using Equation (51):

$$R_{I_n} = 6\frac{D_0}{d^2}n\exp\left(-\frac{E_{\text{diss}}(n)}{kT}\right) \tag{54}$$

## Ostwald ripening



"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
Handcalcs = "e8a07092-c156-4455-ab8e-ed8bc81edefb"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
Match = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
PhysicalConstants = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
CairoMakie = "~0.15.8"
GeometryBasics = "~0.5.10"
Handcalcs = "~0.5.3"
Latexify = "~0.16.10"
Makie = "~0.24.8"
Match = "~2.4.1"
PhysicalConstants = "~0.2.4"
PlutoUI = "~0.7.78"
Symbolics = "~7.0.1"
Unitful = "~1.28.0"
WGLMakie = "~0.13.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "4ceb797e99b5c93a15f7b2d8e5cdb7be37c06cc1"

[[deps.ADTypes]]
git-tree-sha1 = "27cecae79e5cc9935255f90c53bb831cc3c870d7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.18.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bonito]]
deps = ["Base64", "CodecZlib", "Colors", "Dates", "Deno_jll", "HTTP", "Hyperscript", "JSON", "LinearAlgebra", "Markdown", "MbedTLS", "MsgPack", "Observables", "OrderedCollections", "Random", "RelocatableFolders", "SHA", "Sockets", "Tables", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "bb43f72801f703ad3c66833bd02b8f54c7328238"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
version = "4.2.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "5017d6849aff775febd36049f7d926a5fb6677ec"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.8"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "980f01d6d3283b3dbdfd7ed89405f96b7256ad57"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "76dab592fa553e378f9dd8adea16fe2591aa3daa"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.6"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cd6756e833c377e0ce9cd63fb97689a255f12323"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "1.33.4+0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"
weakdeps = ["Makie", "Random"]

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "7bebc8aad6ee6217c78c5ddcf7ed289d65d0263e"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.6"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6522cfb3b8fe97bec632252263057996cbd3de20"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.18.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.Handcalcs]]
deps = ["AbstractTrees", "CodeTracking", "InteractiveUtils", "LaTeXStrings", "Latexify", "MacroTools", "PrecompileTools", "Revise", "TestHandcalcFunctions"]
git-tree-sha1 = "8bcb7f86edf7c6e815667721ea35dbb546ac9dd9"
uuid = "e8a07092-c156-4455-ab8e-ed8bc81edefb"
version = "0.5.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "02b61501dbe6da3b927cc25dacd7ce32390ee970"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.2"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "277779adfedf4a30d66b64edc75dc6bb6d52a16e"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.6"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "4260cfc991b8885bf747801fb60dd4503250e478"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.11"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "e24491cb83551e44a69b9106c50666dea9d953ab"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.4"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "d1b974f376c24dad02c873e951c5cd4e351cd7c2"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.8"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Match]]
deps = ["MacroTools", "OrderedCollections"]
git-tree-sha1 = "58c5c5db26f2f0512facb359991410b7b5982c38"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "2.4.1"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "7eb8cdaa6f0e8081616367c10b31b9d9b34bb02a"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.7"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf"]
git-tree-sha1 = "030f041d5502dbfa41f26f542aaac32bcbe89a64"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.14.0"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsMakieExt = "Makie"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"
weakdeps = ["ChainRulesCore"]

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f2b3b9e52a5eb6a3434c8cca67ad2dde011194f4"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.30+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "df9b7c88c2e7a2e77146223c526bf9e236d5f450"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.PhysicalConstants]]
deps = ["Measurements", "Roots", "Unitful"]
git-tree-sha1 = "9be04ab1a9e54508348a2aacb5b5a4f04c85d397"
uuid = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
version = "0.2.4"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "6122f9423393a2294e26a4efdf44960c5f8acb70"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.78"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "b7e5b731326a99431517b0b4c1f3902e842103a2"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.12.0"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "2f609ec2295c452685d3142bc4df202686e555d2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.16"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5acc6a41b3082920f79ca3c759acbcecf18a8d78"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eee1b9ad8b29ef0d936e3ec9838c7ec089620308"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.16"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "28145feabf717c5d65c1d5e09747ee7b1ff3ed13"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.3"
weakdeps = ["Measurements", "Tables"]

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "49201c2793ce02f141c6f8b5194ce34e8012cd68"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.4"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "30cb5145192c723dff2d5790ca79082f3490079e"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.5.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsChainRulesCoreExt = "ChainRulesCore"
    SymbolicUtilsDistributionsExt = "Distributions"
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "5ccee3582d344a87918840862c0b67285ec9fce1"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.0.1"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestHandcalcFunctions]]
git-tree-sha1 = "54dac4d0a0cd2fc20ceb72e0635ee3c74b24b840"
uuid = "6ba57fb7-81df-4b24-8e8e-a3885b6fcae7"
version = "0.2.4"

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WGLMakie]]
deps = ["Bonito", "Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "LinearAlgebra", "Makie", "Observables", "PNGFiles", "PrecompileTools", "RelocatableFolders", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "32801246eb6c7afb0e1e49509b3ffebecb538657"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.13.8"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e015f211ebb898c8180887012b938f3851e719ac"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.55+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╟─f5bc5f70-0a80-11f1-a259-8be415f4d105
# ╟─6d0c96f2-4dfd-43a1-86af-e6aaa9b04f0c
# ╟─4b5af7d2-3894-4220-8f2d-dec15268305d
# ╟─0dc3df3c-0999-4f9c-8e68-4b71a51d00ed
# ╟─45f0efff-bc20-49d3-8442-3948d82b6441
# ╟─8b4c2f8f-6027-43f9-8bc2-d2fbd7cbf3ea
# ╟─c8a756c6-14d3-4bb9-ab4f-e794a0653344
# ╟─1be863ee-0514-4de1-a0c4-9460cf2852c9
# ╟─a0a819cd-7d17-4d61-a88d-a2a9fdd3e8ab
# ╟─48a4e303-921c-4b59-93a7-633f2f2cffb8
# ╟─14b9e9d3-a71b-42a4-9942-94ebfe63ccb4
# ╟─db101921-c3dd-4d69-956a-5bc6650dd273
# ╟─3e1815e7-8306-4aa3-ada8-f113d2f4785f
# ╟─670b6029-7c3d-4f9d-a218-1dbf402667d0
# ╟─626913f7-4038-4af0-bf7d-974a85f9cffc
# ╟─d23e4178-1741-4619-a54c-7d505996ff5e
# ╟─d1766538-b8d7-458b-951e-01b23a3924cd
# ╟─51869391-d0e7-408a-9d8d-873afa01976c
# ╟─3428b15a-6d55-4bcf-b361-00af9d3de667
# ╟─8518ab0e-173f-4ee2-a0a8-7e3eeaf4fcc6
# ╟─ba5df24e-436b-4269-b355-8d52e7a00be2
# ╟─63cd3890-db15-4239-9fae-0e8565f8a7ac
# ╟─891ace0a-83b4-4f85-b571-ad19969aea22
# ╟─c1a1aa76-1b8d-4abf-a697-83d9ffab3875
# ╟─8f6d9c95-0cfa-41c1-ae3c-c5e7ec6d8fd4
# ╟─afc0732f-f49a-4b27-ade6-88e0936ccd5a
# ╟─ab388d5e-651d-4962-9b7e-1dfdf39b35fe
# ╟─7c249d66-6ab2-4161-977e-abc9f719073a
# ╟─b0c71984-34cc-4496-a1a9-e2a01ac1dea7
# ╟─0188a13e-0b79-496c-992b-74227047529f
# ╟─dcd81bf5-c7e4-4c93-a7f4-6657d638e0b4
# ╟─34c17de5-576a-40ba-b5a4-afa9ee2d8aef
# ╟─5b23e294-6bf6-439c-84ff-34cf99d31e53
# ╟─c1feebf7-4384-4196-bede-8f0662f301c4
# ╟─9ae89154-24de-42b6-81ae-ebb07d6a87ea
# ╟─9ca2fdcf-5909-4c52-b470-080add85212a
# ╟─dc173930-d763-4e53-baca-be8779f41614
# ╟─591e64a6-56d6-4014-8db6-d90274f5216f
# ╟─2322c0d0-3fa4-47fd-8c3c-905af02461d5
# ╟─16d9d932-ddf4-4486-8324-ae6dcdf86e18
# ╟─4365438f-afa6-4478-9d86-6c82545dc84e
# ╟─f5ef096e-de1f-4299-a8d0-4f9d757f96a2
# ╟─a1b55842-ef12-44f0-a686-4976cd55765a
# ╟─7b0f6bf6-eb66-4f1a-bc5c-71cfbe2b0c90
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
