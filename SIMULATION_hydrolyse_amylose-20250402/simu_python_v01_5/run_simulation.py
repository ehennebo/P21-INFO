# Framework for simulation of molecular diffusion - INSA
# Catherine Pothier and Christophe Rigotti

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

from Simulator import Simulator


nbIterations = 200
simu = Simulator(nbEnzymes=100, nbPolymers=100,
             nbIterations=nbIterations, xMax=80, yMax=80)

time_start = time.time()
print("Start simulation")
universeViews = simu.run()
print(f"Simulation duration: {time.time()-time_start}")


print("Start computation of animation")
# Building an animation using the universe views
firstFigure = plt.figure("Simulation")
axes = firstFigure.add_subplot()

allSnapshots = []
plt.axis("off")

for i in range(nbIterations):    
    snapshot = axes.imshow(universeViews[i,:,:], animated=True)
    allSnapshots.append([snapshot])

anim = animation.ArtistAnimation(firstFigure, allSnapshots, interval=20,
                                    blit=True, repeat=False)

print("End of computation of animation")
# To display the animation (comment the line if you don't want to see it)
plt.show()
print("plt.show done")
# CAUTION: depending on the execution mode (IDE, console, ...),
# may required to close the figure window to go on ...
# in some mode the animation can be shown only at the end
# of the execution

save_animation = False
if save_animation:
    # To write the animation to a gif file:
    print("Start writing of animation in file")
    gifWriter = animation.PillowWriter(fps=30)
    anim.save("out_animation.gif", writer=gifWriter)
    print("End of animation writing")


print("Bye ...")
