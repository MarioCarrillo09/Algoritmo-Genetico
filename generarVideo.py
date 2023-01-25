import cv2
import numpy as np

# choose codec according to format needed
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter('video.avi', fourcc, 1, (640, 480))

for j in range(1, 21):
    img = cv2.imread("./images/generation " + str(j) + '.png')
    video.write(img)

cv2.destroyAllWindows()
video.release()
