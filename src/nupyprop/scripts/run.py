#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 21:02:17 2020

@author: sam
"""

import argparse
import nupyprop.main as Main

import time
import numpy as np


def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-e",
        "--energy",
        dest="energy_val",
        nargs="?",
        const=f"{np.array2string(np.linspace(7, 11, 17), separator=',')}"[1:-1],
        default=f"{np.array2string(np.linspace(7, 11, 17), separator=',')}"[1:-1],
        help="log_10 value of incoming neutrino energy; defaults are 7-11 GeV in quarter"
        " decades",
    )

    parser.add_argument(
        "-a",
        "--angle",
        dest="angle_val",
        nargs="?",
        const=f"{np.array2string(np.arange(1, 36), separator=',')}"[1:-1],
        default=f"{np.array2string(np.arange(1, 36), separator=',')}"[1:-1],
        help="value of Earth emergence angle; defaults are 1-35 degrees, in steps of 1 degree.",
    )

    parser.add_argument(
        "-i",
        "--idepth",
        dest="idepth_val",
        nargs="?",
        type=float,
        const=4,
        default=4,
        help="value of water layer in km; default is 4 km",
    )

    parser.add_argument(
        "-l",
        "--lepton",
        dest="lepton_id",
        nargs="?",
        type=str,
        const="tau",
        default="tau",
        help="particle for energy loss and propagation - can be tau or muon; default is tau",
    )

    parser.add_argument(
        "-n",
        "--nu_type",
        dest="nu_type_id",
        nargs="?",
        type=str,
        const="neutrino",
        default="neutrino",
        help="type of neutrino matter - can be neutrino or anti_neutrino; default is neutrino",
    )

    parser.add_argument(
        "-t",
        "--energy_loss",
        dest="loss_type",
        nargs="?",
        type=str,
        const="stochastic",
        default="stochastic",
        help="energy loss type for lepton - can be stochastic or continuous; default is stochastic",
    )

    parser.add_argument(
        "-x",
        "--xc_model",
        dest="xc_model_id",
        nargs="?",
        type=str,
        const="ct18nlo",
        default="ct18nlo",
        help="neutrino cross-section model; default is ct18nlo",
    )

    parser.add_argument(
        "-p",
        "--pn_model",
        dest="pn_model_id",
        nargs="?",
        type=str,
        const="allm",
        default="allm",
        help="lepton photonuclear energy loss model; default is allm",
    )

    parser.add_argument(
        "-f",
        "--fac_nu",
        dest="fac_nu_val",
        nargs="?",
        type=float,
        const=1.0,
        default=1.0,
        help="rescaling for BSM neutrino cross-sections; default is 1.0",
    )

    parser.add_argument(
        "-s",
        "--stat",
        dest="stat_val",
        nargs="?",
        type=float,
        const=1e7,
        default=1e7,
        help="statistics or number of neutrinos injected; default is 1e7",
    )

    parser.add_argument(
        "-c",
        "--cdf_only",
        dest="cdf_id",
        nargs="?",
        type=str,
        const="no",
        default="no",
        help="CDF only; default is no. If set to yes, the output file will NOT contain outgoing lepton energies.",
    )

    return parser


def main():

    parser = get_parser()
    args = parser.parse_args()

    energies = np.power(10, np.fromstring(args.energy_val, dtype=float, sep=","))
    angles = np.fromstring(args.angle_val, dtype=int, sep=",")

    idepth = int(args.idepth_val)
    lepton = str(args.lepton_id)
    nu_type = str(args.nu_type_id)
    type_loss = str(args.loss_type)
    cross_section_model = str(args.xc_model_id)
    pn_model = str(args.pn_model_id)
    fac_nu = float(args.fac_nu_val)
    stat = int(args.stat_val)
    cdf_only = str(args.cdf_id)

    start_time = time.time()

    Main.main(
        energies,
        angles,
        nu_type,
        cross_section_model,
        pn_model,
        idepth,
        lepton,
        fac_nu,
        stat,
        type_loss,
        cdf_only,
    )

    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute")
