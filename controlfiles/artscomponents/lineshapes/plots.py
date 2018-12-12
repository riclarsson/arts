#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:14:03 2018

Plots testdata/dtest-"+ls+"/ files.

@author: larsson
"""

import typhon
import matplotlib.pyplot as plt
import numpy as np

""" NOTE:  fake-htp gives bad VMR comparison before internal code is fixed """

ls = "lorentz"
pm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat.xml")
adpm = typhon.arts.xml.load("testdata/test-"+ls+"/dpropmat.xml")
f = np.linspace(90, 110, 1001)
plt.figure(figsize=(5, 5))
plt.plot(f, pm[0].data[0, 0, :, 0])
plt.show()

plt.figure(figsize=(15, 15))
i = 1
for ipm in adpm:
    plt.subplot(6, 6, i)
    plt.plot(f, abs(ipm.data[0, 0, :, 0]))
    i += 1
plt.tight_layout()

pmd =typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dT.xml")
plt.subplot(6, 6, 1)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.0001)
plt.title("Temperature")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-df.xml")
plt.subplot(6, 6, 2)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 100)
plt.title("Frequency")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dvmr.xml")
plt.subplot(6, 6, 3)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.0001)
plt.title("VMR")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-ds0.xml")
plt.subplot(6, 6, 4)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-30)
plt.title("Line Strength")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-df0.xml")
plt.subplot(6, 6, 5)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 10)
plt.title("Line Center")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-G0-X0.xml")
plt.subplot(6, 6, 6)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 10)
plt.title("SELF-G0-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-G0-X1.xml")
plt.subplot(6, 6, 7)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("SELF-G0-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-D0-X0.xml")
plt.subplot(6, 6, 9)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 10)
plt.title("SELF-D0-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-D0-X1.xml")
plt.subplot(6, 6, 10)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("SELF-D0-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-Y-X0.xml")
plt.subplot(6, 6, 12)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-10)
plt.title("SELF-Y-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-Y-X1.xml")
plt.subplot(6, 6, 13)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-12)
plt.title("SELF-Y-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-Y-X2.xml")
plt.subplot(6, 6, 14)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("SELF-Y-X2")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-G-X0.xml")
plt.subplot(6, 6, 15)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-14)
plt.title("SELF-G-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-G-X1.xml")
plt.subplot(6, 6, 16)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-16)
plt.title("SELF-G-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-G-X2.xml")
plt.subplot(6, 6, 17)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("SELF-G-X2")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-DV-X0.xml")
plt.subplot(6, 6, 18)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.1)
plt.title("SELF-DV-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-DV-X1.xml")
plt.subplot(6, 6, 19)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.01)
plt.title("SELF-DV-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-SELF-DV-X2.xml")
plt.subplot(6, 6, 20)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("SELF-DV-X2")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-G0-X0.xml")
plt.subplot(6, 6, 21)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 10)
plt.title("AIR-G0-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-G0-X1.xml")
plt.subplot(6, 6, 22)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("AIR-G0-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-D0-X0.xml")
plt.subplot(6, 6, 24)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 10)
plt.title("AIR-D0-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-D0-X1.xml")
plt.subplot(6, 6, 25)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("AIR-D0-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-Y-X0.xml")
plt.subplot(6, 6, 27)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-10)
plt.title("AIR-Y-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-Y-X1.xml")
plt.subplot(6, 6, 28)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-12)
plt.title("AIR-Y-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-Y-X2.xml")
plt.subplot(6, 6, 29)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("AIR-Y-X2")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-G-X0.xml")
plt.subplot(6, 6, 30)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-14)
plt.title("AIR-G-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-G-X1.xml")
plt.subplot(6, 6, 31)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-16)
plt.title("AIR-G-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-G-X2.xml")
plt.subplot(6, 6, 32)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("AIR-G-X2")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-DV-X0.xml")
plt.subplot(6, 6, 33)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.1)
plt.title("AIR-DV-X0")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-DV-X1.xml")
plt.subplot(6, 6, 34)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.01)
plt.title("AIR-DV-X1")

pmd = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlf-AIR-DV-X2.xml")
plt.subplot(6, 6, 35)
plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-05)
plt.title("AIR-DV-X2")

plt.tight_layout()
plt.show()