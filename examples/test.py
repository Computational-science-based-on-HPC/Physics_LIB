import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys
# Constants
ball_radius = 0.1  # Radius of the ball
anchor_x = 0.0  # x-coordinate of the anchor point
anchor_y = 10.0  # y-coordinate of the anchor point

# Read relative position file
# Replace with the actual path to the position file
input_file = "/media/slinga/HDD/FCAI/GP/Codes/LIB/lib/Physics_LIB/examples/Simulation/Damped oscillation serial/damped_os_serial_displacement_2023-07-13_23-42-04.sim"

relative_positions = np.loadtxt(input_file)
frames = 16*24
step_size = len(relative_positions)//frames
relative_positions = relative_positions[::step_size]
# Function to update the frame


def update_frame(frame):
    print(frame)
    plt.cla()  # Clear the current plot

    # Calculate ball position
    ball_x = anchor_x
    ball_y = anchor_y-relative_positions[frame]

    # Plot the anchor point
    plt.plot(anchor_x, anchor_y, 'ko')

    # Plot the ball
    plt.plot(ball_x, ball_y, 'ro', markersize=ball_radius*100)
    plt.plot([anchor_x, ball_x], [anchor_y, ball_y], 'k--')
    # Set plot limits
    plt.xlim(anchor_x - 2*ball_radius, anchor_x + 2*ball_radius)
    plt.ylim(max(relative_positions) - 2*ball_radius,
             anchor_y-min(relative_positions) + 2*ball_radius)

    # Add labels and title
    plt.xlabel('Position')
    plt.ylabel('Height')
    plt.title('Damped Oscillation')
    if frame == len(relative_positions) - 1:
        sys.exit()


# Create animation
fig = plt.figure()

animation = animation.FuncAnimation(fig, update_frame,
                                    frames=len(relative_positions),
                                    interval=42)

# # Save animation as a GIF file
# Replace with the desired output path and filename
output_file = "/media/slinga/HDD/FCAI/GP/Codes/LIB/lib/Physics_LIB/graphs/animation.gif"
animation.save(output_file, writer='pillow')


# Save animation as a GIF file using 'imagemagick' writer
# output_file = "/media/slinga/HDD/FCAI/GP/Codes/LIB/lib/Physics_LIB/graphs/animation.gif"  # Replace with the desired output path and filename
# animation.save(output_file, writer='imagemagick')
