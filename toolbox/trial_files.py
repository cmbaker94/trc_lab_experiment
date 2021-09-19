# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 12:41:18 2021

@author: cmbaker9

This function stores all the trial information for all processed cases.
"""
def trial_files(Hs,Tp,spread,h):
    if spread == 0:
        if Hs == 0.3:
            Tinfo = {
                "Hs": 0.3,
                "Tp": 2,
                "h": 1.07,
                "spread": 0,
                "day": '7',
                "trial": '14',
                "insitu": "09-01-2018-2213UTC",
                "clpath": "09-01-2018-2214UTC",
                "timesection": "2229-2238",
                "camera": "09-01-2018-2155UTC_Scene1",
                "frames": "07200-11999",
                "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                "parray":'sz'
                }
        elif Hs == 0.25:
            Tinfo = {
                "Hs": 0.25,
                "Tp": 2,
                "h": 1.07,
                "spread": 0,
                "day": '7',
                "trial": '12',
                "insitu": "09-01-2018-2001UTC",
                "clpath": "09-01-2018-2002UTC",
                "timesection": "2017-2026",
                "camera": "09-01-2018-1950UTC_Scene1",
                "frames": "07200-11999",
                "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                "parray":'sz'
                }
    elif spread == 10:
        Tinfo = {
                "Hs": 0.25,
                "Tp": 2,
                "h": 1.07,
                "spread": 10,
                "day": '5',
                "trial": '05',
                "insitu": "08-30-2018-2026UTC",
                "clpath": "09-06-2018-2342UTC",
                "timesection": "2357-0006",
                "camera": "09-06-2018-2334UTC_Scene1",
                "frames": "07200-11999",
                "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                "parray":'sz'
                }
    elif spread == 20:
        if Hs == 0.3:
            if Tp == 2:
                 Tinfo = {
                 "Hs": 0.3,
                 "Tp": 2,
                 "h": 1.07,
                 "spread": 20,
                 "day": '5',
                 "trial": '07',
                 "insitu": "08-30-2018-2222UTC",
                 "clpath": "09-06-2018-1655UTC",
                 "timesection": "1710-1719",
                 "camera": "09-06-2018-1646UTC_Scene1",
                 "frames": "07200-11999",
                 "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                 "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                 "parray":'sz'
                 }
            elif Tp == 3:
                Tinfo = {
                    "Hs": 0.3,
                    "Tp": 3,
                    "h": 1.07,
                    "spread": 20,
                    "day": '7',
                    "trial": '09',
                    "insitu": "09-01-2018-1626UTC",
                    "clpath": "09-01-2018-1627UTC",
                    "timesection": "1642-1651",
                    "camera": "09-01-2018-1530UTC_Scene1",
                    "frames": "07200-11999",
                    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                    "parray":'sz'
                    }
        elif Hs == 0.25:
             Tinfo = {
                    "Hs": 0.25,
                    "Tp": 2,
                    "h": 1.07,
                    "spread": 20,
                    "day": '5',
                    "trial": '04',
                    "insitu": "08-30-2018-1904UTC",
                    "clpath": "09-06-2018-2049UTC",
                    "timesection": "2104-2113",
                    "camera": "09-06-2018-2043UTC_Scene1",
                    "frames": "07200-11999",
                    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                    "parray":'sz'
                    }
    elif spread == 30:
        if Hs == 0.3:
            Tinfo = {
            "Hs": 0.3,
            "Tp": 2,
            "h": 1.07,
            "spread": 30,
            "day": '4',
            "trial": '06',
            "insitu": "08-29-2018-2255UTC",
            "clpath": "08-29-2018-2255UTC",
            "timesection": "2310-2319",
            "camera": "08-29-2018-2236UTC_Scene1",
            "frames": "07200-11999",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'sz'
            }
        elif Hs == 0.25:
            Tinfo = {
            "Hs": 0.25,
            "Tp": 2,
            "h": 1.07,
            "spread": 30,
            "day": '7',
            "trial": '13',
            "insitu": "09-01-2018-2102UTC",
            "clpath": "09-01-2018-2102UTC",
            "timesection": "2117-2126",
            "camera": "09-01-2018-2052UTC_Scene1",
            "frames": "07200-11999",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'sz'
            }
    elif spread == 40:
        if Hs == 0.3:
            if h == 1.07:
                if Tp == 2:
                    Tinfo = {
                        "Hs": 0.3,
                        "Tp": 2,
                        "h": 1.07,
                        "spread": 40,
                        "day": '4',
                        "trial": '07',
                        "insitu": "08-29-2018-2358UTC",
                        "clpath": "08-29-2018-2359UTC",
                        "timesection": "0014-0023",
                        "camera": "08-29-2018-2349UTC_Scene1",
                        "frames": "07200-11999",
                        "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                        "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                        "parray":'sz'
                        }
                elif Tp == 3:
                    Tinfo = {
                        "Hs": 0.3,
                        "Tp": 3,
                        "h": 1.07,
                        "spread": 40,
                        "day": '7',
                        "trial": '10',
                        "insitu": "09-01-2018-1757UTC",
                        "clpath": "09-01-2018-1758UTC",
                        "timesection": "1813-1822",
                        "camera": "09-01-2018-1740UTC_Scene1",
                        "frames": "07200-11999",
                        "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                        "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                        "parray":'sz'
                        }
            elif h == 1.00:
                Tinfo = {
                    "Hs": 0.3,
                    "Tp": 2,
                    "h": 1.00,
                    "spread": 40,
                    "day": '6',
                    "trial": '07',
                    "insitu": "08-31-2018-2232UTC",
                    "clpath": "08-31-2018-2232UTC",
                    "timesection": "2247-2256",
                    "camera": "08-31-2018-2225UTC_Scene1",
                    "frames": "07200-11999",
                    "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
                    "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                    "parray":'sz'
                    }
        elif Hs == 0.25:
            Tinfo = {
            "Hs": 0.25,
            "Tp": 2,
            "h": 1.07,
            "spread": 40,
            "day": '5',
            "trial": '02',
            "insitu": "08-30-2018-1655UTC",
            "clpath": "08-30-2018-1655UTC",
            "timesection": "1710-1719",
            "camera": "08-30-2018-1634UTC_Scene1",
            "frames": "07200-11999",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'sz'
            }
        elif Hs == 0.2:
            Tinfo = {
            "Hs": 0.2,
            "Tp": 2,
            "h": 1.07,
            "spread": 40,
            "day": '5',
            "trial": '01',
            "insitu": "08-30-2018-1534UTC",
            "clpath": "08-30-2018-1534UTC",
            "timesection": "1549-1558",
            "camera": "08-30-2018-1518UTC_Scene1",
            "frames": "07200-11999",
            "wg": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
            "press":[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            "parray":'sz'
            }
    return Tinfo