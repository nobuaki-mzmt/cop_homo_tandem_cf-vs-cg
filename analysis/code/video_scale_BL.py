"""
Program to measure scale / body length from video files manually
videoScaleBL.py
N. Mizumoto
"""

# ---------------------------------------------------------------------------------------------
# packages
# ---------------------------------------------------------------------------------------------
import PySimpleGUI as sg
import cv2
import sys
import numpy as np
import glob
import tqdm
from keyboard import press
import os
import math
import pandas as pd
import pickle
from numpy.linalg import norm


# ---------------------------------------------------------------------------------------------
# Main part of video analysis
# ---------------------------------------------------------------------------------------------
def ImageAnalysis(idir, odir, img_scale, measure_scale, shape, num_ind, frame_interval, skip_analyzed):

    # Loading data
    path = glob.glob(idir + os.sep + "*.mp4")
    file_nums = list(range((len(path))))

    # Dataframe
    if os.path.exists(odir + "/res.pickle"):
        print("existing analysis loaded")
        with open(odir+ os.sep  + 'res.pickle', mode='rb') as f:
            df = pickle.load(f)
    else:
        df_column = ["name", "width", "height", "length", "fps", "frame", "scale"]
        for i in range(num_ind):
            df_column.append("bodyLength" + str(i))
        df = pd.DataFrame(np.arange(len(path)*(7+num_ind)). reshape(len(path), 7+num_ind),
                          columns=df_column)

    # ------------------ main ------------------
    def frame_check(event, x, y, flags, param):
            if event == cv2.EVENT_LBUTTONDOWN:
                press('l')
            if event == cv2.EVENT_RBUTTONDOWN:
                press('r')

    for i in file_nums:
        cv2.namedWindow(winname='window')

        # ----- file info -----
        v = path[i]
        file_name = os.path.basename(v)
        name = file_name.split('.')[0]
        if df.iloc[i, 0] == name:
            if skip_analyzed == "true":
                continue
            img_read = cv2.imread(odir + name + ".jpg")
            cv2.putText(img_read, 'Skip? R:yes, L:No',
                        (10, 100), cv2.FONT_HERSHEY_PLAIN, 2, (0, 0, 255), 2, cv2.LINE_AA)
            cv2.imshow('window', img_read)
            cv2.setMouseCallback('window', frame_check)
            k = cv2.waitKey(0)
            if k == 27:
                break
            if k == ord("r"):
                continue

        video = cv2.VideoCapture(v)
        width = video.get(cv2.CAP_PROP_FRAME_WIDTH)
        height = video.get(cv2.CAP_PROP_FRAME_HEIGHT)
        count = video.get(cv2.CAP_PROP_FRAME_COUNT)
        fps = video.get(cv2.CAP_PROP_FPS)
        print("name:{}, width:{}, height:{}, count:{}, fps:{}".format(name, width, height, count, fps))
        img_shape = np.array([int(width), int(height)])
        # ----------

        # region ----- 1. Extract frame -----
        img = None
        frame_id = 0

        while True:
            video.set(cv2.CAP_PROP_POS_FRAMES, frame_id)
            ret, frame = video.read()
            frame = cv2.resize(frame, dsize=None, fx=img_scale, fy=img_scale)
            frame_copy = frame.copy()
            cv2.putText(frame_copy, 'file: ' + name + '  frame: ' + str(frame_id),
                        (10, 50), cv2.FONT_HERSHEY_PLAIN, 2, (0, 0, 255), 2, cv2.LINE_AA)
            cv2.putText(frame_copy, 'Use this frame? R:yes, L:No',
                        (10, 100), cv2.FONT_HERSHEY_PLAIN, 2, (0, 0, 255), 2, cv2.LINE_AA)
            cv2.imshow('window', frame_copy)
            cv2.setMouseCallback('window', frame_check)
            k = cv2.waitKey(0)
            if k == ord("r"):
                img = frame.copy()
                break
            elif k == ord("l"):
                frame_id   = frame_id + frame_interval
                if frame_id > count:
                    print("End of frames. Maybe reduce frame_interval to sample more frames.")
                    frame_id = 0
            elif k == 27:
                break

        if k == 27:
            break
        # endregion ----------

        img_output = img.copy()

        # region ----- 2. Scale -----
        if measure_scale == "True":
            sx0, sy0 = 0, 0
            sx1, sy1 = 0, 0
            drawing = False

            if shape == "circle":
                def scale_draw(event, x, y, flags, param):
                    nonlocal sx0, sy0, scale, drawing, end, img_output, img_copy
                    if event == cv2.EVENT_LBUTTONDOWN:
                        drawing = True
                        sx0, sy0 = x, y
                    elif event == cv2.EVENT_MOUSEMOVE:
                        if drawing:
                            cv2.circle(img_copy, center=(int((x - sx0) / 2) + sx0, int((y - sy0) / 2) + sy0),
                                       radius=int(math.sqrt((sx0 - x) ** 2 + (sy0 - y) ** 2) / 2),
                                       color=(0, 0, 255), thickness=1)
                    elif event == cv2.EVENT_LBUTTONUP:
                        scale = math.sqrt((sx0 - x) ** 2 + (sy0 - y) ** 2)
                        cv2.circle(img_copy, center=(int((x - sx0) / 2) + sx0, int((y - sy0) / 2) + sy0),
                                   radius=int(scale / 2),
                                   color=(0, 0, 255), thickness=1)
                        cv2.circle(img_output, center=(int((x - sx0) / 2) + sx0, int((y - sy0) / 2) + sy0),
                                   radius=int(scale / 2),
                                   color=(0, 0, 255), thickness=1)

                        drawing = False
                    if event == cv2.EVENT_RBUTTONDOWN:
                        end = 1
                    press('enter')

            elif shape == "line":
                def scale_draw(event, x, y, flags, param):
                    nonlocal sx0, sy0, scale, drawing, end, img_output
                    if event == cv2.EVENT_LBUTTONDOWN:
                        drawing = True
                        sx0, sy0 = x, y
                    elif event == cv2.EVENT_MOUSEMOVE:
                        if drawing:
                            cv2.line(img_copy, (sx0, sy0), (x, y), (0, 0, 255), 1)
                    elif event == cv2.EVENT_LBUTTONUP:
                        cv2.line(img_copy, (sx0, sy0), (x, y), (0, 0, 255), 1)
                        cv2.line(img_output, (sx0, sy0), (x, y), (0, 0, 255), 1)
                        drawing = False
                        scale = math.sqrt((sx0 - x) ** 2 + (sy0 - y) ** 2)
                    if event == cv2.EVENT_RBUTTONDOWN:
                        end = 1


            img_copy = img.copy()
            end = 0
            drawing = False
            cv2.putText(img_copy, 'Scaling', (10, 50), cv2.FONT_HERSHEY_PLAIN, 1, (0, 0, 255), 1, cv2.LINE_AA)
            while True:
                cv2.imshow('window', img_copy)
                if drawing:
                    img_copy = img.copy()
                cv2.setMouseCallback('window', scale_draw)
                k = cv2.waitKey(0)
                if end == 1:
                    break
                if k == 27:
                    break
            if k == 27:
                break
        
        else:
            scale = 0
        # endregion

        drawing = False

        # region --- 3.  Measure body length --- #
        img_copy = img.copy()
        cv2.putText(img_copy, 'Body Length', (10, 50), cv2.FONT_HERSHEY_PLAIN, 1, (0, 0, 255), 1, cv2.LINE_AA)
        body_length = [0] * num_ind
        bl_data = []

        # todo: code for zooming is not great. but I have no idea how to improve yet.
        zoom, zoom_xy = [1, 1, 1], [np.array([0, 0]), np.array([0, 0]), np.array([0, 0])]  # for x2, x4, x8
        mouse_xy = np.array([0, 0])

        for ii in range(num_ind):
            current_bl = [[0, 0], [0, 0]]
            bl = 0
            while True:
                cv2.imshow('window', img_copy)
                def bl_draw(event, x, y, flags, param):
                    nonlocal bl_data, ii, img, img_copy, current_bl, mouse_xy, zoom, zoom_xy, bl

                    if event == cv2.EVENT_MOUSEMOVE:
                        mouse_xy = np.array([x, y])

                    if event == cv2.EVENT_LBUTTONDOWN:
                        current_bl[0] = (((np.array([x, y])+zoom_xy[2])/zoom[2]+zoom_xy[1])/zoom[1]+zoom_xy[0])/zoom[0]
                        current_bl[0] = current_bl[0].astype(int)

                    if event == cv2.EVENT_LBUTTONUP:
                        current_bl[1] = (((np.array([x, y])+zoom_xy[2])/zoom[2]+zoom_xy[1])/zoom[1]+zoom_xy[0])/zoom[0]
                        current_bl[1] = current_bl[1].astype(int)
                        cv2.line(img_copy, 
                            ((current_bl[0]*zoom[0]-zoom_xy[0])*zoom[1]-zoom_xy[1])*zoom[2]-zoom_xy[2],
                            ((current_bl[1]*zoom[0]-zoom_xy[0])*zoom[1]-zoom_xy[1])*zoom[2]-zoom_xy[2], (0, 255, 255), 2)
                        bl = norm(current_bl[1] - current_bl[0])
                    if event == cv2.EVENT_RBUTTONDOWN:
                        cv2.line(img_copy, 
                            ((current_bl[0]*zoom[0]-zoom_xy[0])*zoom[1]-zoom_xy[1])*zoom[2]-zoom_xy[2],
                            ((current_bl[1]*zoom[0]-zoom_xy[0])*zoom[1]-zoom_xy[1])*zoom[2]-zoom_xy[2], (0, 0, 255), 2)
                        cv2.line(img_output, current_bl[0], current_bl[1], (0, 0, 255), 2)
                        
                        press("r")
                    cv2.imshow('window', img_copy)

                cv2.setMouseCallback('window', bl_draw)
                k = cv2.waitKey(0)

                def zoom_func(img_z, mouse_xy, img_shape, zoom):
                    mouse_xy[0] = max(mouse_xy[0], img_shape[0] / 4)
                    mouse_xy[1] = max(mouse_xy[1], img_shape[1] / 4)
                    mouse_xy[0] = min(mouse_xy[0], img_shape[0] * 3 / 4)
                    mouse_xy[1] = min(mouse_xy[1], img_shape[1] * 3 / 4)
                    img_zoom = cv2.resize(img_z, dsize=(img_shape * 2))
                    img_zoom = img_zoom[int(mouse_xy[1] * 2 - img_shape[1] / 2):int(mouse_xy[1] * 2 + img_shape[1] / 2),
                               int(mouse_xy[0] * 2 - img_shape[0] / 2):int(mouse_xy[0] * 2 + img_shape[0] / 2)]
                    zoom_xy = mouse_xy * 2 - img_shape / 2
                    zoom_xy = zoom_xy.astype(int)
                    return img_zoom, zoom_xy, zoom*2

                if k == ord("z"):
                    if zoom[1] == 2 and zoom[2] == 1:
                        img_copy, zoom_xy[2], zoom[2] = zoom_func(img_copy, mouse_xy, img_shape, zoom[2])
                    elif zoom[0] == 2 and zoom[1] == 1:
                        img_copy, zoom_xy[1], zoom[1] = zoom_func(img_copy, mouse_xy, img_shape, zoom[1])
                    elif zoom[0] == 1:
                        img_copy, zoom_xy[0], zoom[0] = zoom_func(img_copy, mouse_xy, img_shape, zoom[0])
                
                if k == ord("x"):
                    # cancel zoom when redo
                    img_copy = img.copy()
                    if ii > 0:
                        for i in range(ii-1):
                            cv2.line(img_copy, 
                                ((bl_data[ii][0]*zoom[0]-zoom_xy[0])*zoom[1]-zoom_xy[1])*zoom[2]-zoom_xy[2],
                                ((bl_data[ii][1]*zoom[0]-zoom_xy[0])*zoom[1]-zoom_xy[1])*zoom[2]-zoom_xy[2], (0, 0, 255), 2)
                        
                    zoom, zoom_xy = [1, 1, 1], [np.array([0, 0]), np.array([0, 0]), np.array([0, 0])]

                if k == ord("r"):
                    bl_data.append(current_bl)
                    body_length[ii] = bl
                    break

                if k == 27:
                    break
            if k == 27:
                break
        if k == 27:
            break

        # endregion



        # region ----- 4. Output -----#
        df.iloc[i:(i+1), 0:7] = [name, width, height, count, fps, frame_id, scale/img_scale,]
        for ii in range(num_ind):
            df.iloc[i:(i+1), 7+ii] = body_length[ii]/img_scale
        cv2.putText(img_output, name, (10, 50), cv2.FONT_HERSHEY_PLAIN, 2, (0, 0, 255), 2, cv2.LINE_AA)
        cv2.imwrite(odir + os.sep + name + ".jpg", img_output)

    cv2.destroyAllWindows()

    print(df)

    with open(odir + os.sep + 'res.pickle', mode='wb') as f:
        pickle.dump(df, f)

    df.to_csv(odir + os.sep + '/res.csv')


# ---------------------------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------------------------
sg.theme('Dark')
frame1 = sg.Frame('', [
    [sg.Text("In"),
     sg.InputText('Input folder', key='-INPUTTEXT-', enable_events=True, ),
     sg.FolderBrowse(button_text='select', size=(8, 1), key="-INFOLDERNAME-")
     ],
    [sg.Text("Out"),
     sg.InputText('Output folder', key='-INPUTTEXT-', enable_events=True, ),
     sg.FolderBrowse(button_text='select', size=(8, 1), key="-OUTFOLDERNAME-"),
     sg.Text("(* will be created if not specified)")
     ],
     [sg.Text("Display img relative size (default = 0.5)"),
     sg.In(key='-scale-')],
     [sg.Text("Measure scale"),
     sg.Combo(['True', 'False'], default_value="True", size=(4, 1), key="-measure_scale")
     ],     
     [sg.Text("Methods"),
     sg.Combo(['circle', 'line'], default_value="circle", size=(4, 1), key="-shape")
     ],
     [sg.Text("Num of individual (default = 2)"),
     sg.In(key='-NUMIND-')],
     [sg.Text("Sample frame interval (default = 3000)"),
     sg.In(key='-FRAMEINTERVAL-')],
     [sg.Text("skip analyzed files", size=(15,1)),
              sg.Combo(['true', 'false'], default_value="true", size=(6, 1), key="-SKIP_ANALYZED-")],
], size=(800, 300))

frame2 = sg.Frame('', [
    [sg.Submit(button_text='Start', size=(10, 3), key='button_start')]], size=(100, 140))

layout = [[frame1, frame2]]

window = sg.Window('files', layout, resizable=True)


while True:
    event, values = window.read()

    if event is None:
        print('exit')
        break

    if len(values["-INFOLDERNAME-"]) == 0:
            print("no input!")
    else:
        if event == 'button_start':
            idir = values["-INFOLDERNAME-"]
            if len(values["-OUTFOLDERNAME-"]) == 0:
                odir = idir + "/ouput/"
                if not os.path.exists(odir):
                    os.makedirs(odir)
            else:
                odir = values["-OUTFOLDERNAME-"]

            if len(values["-scale-"]) == 0:
                img_scale = 0.5
            else:
                img_scale = float(values["-scale-"])

            if len(values["-NUMIND-"]) == 0:
                num_ind = 2
            else:
                num_ind = int(values["-NUMIND-"])

            if len(values["-FRAMEINTERVAL-"]) == 0:
                frame_interval = 3000
            else:
                frame_interval = int(values["-FRAMEINTERVAL-"])
            skip_analyzed = values["-SKIP_ANALYZED-"]
            shape = values["-shape"]    
            measure_scale = values["-measure_scale"]    

            ImageAnalysis(idir, odir, img_scale, measure_scale, shape, num_ind, frame_interval, skip_analyzed)
window.close()



