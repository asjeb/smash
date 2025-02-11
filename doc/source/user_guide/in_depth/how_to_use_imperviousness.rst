.. _user_guide.in_depth.how_to_use_imperviousness:

=========================
How To Use Imperviousness
=========================

This tutorial explains how to use a dataset of imperviousness. 
The coefficients of imperviousness are applied as area percentage on each pixel.
The user must supply to the setup a `.tif` file at the same resolution as the SMASH grid.
The user must add the following lines to the setup:

.. code-block:: python

    "read_imperviousness":True,
    "imperviousness_format":"tif",
    "imperviousness_file": "./imperviousness.tif",

We suggest to read a map of imperviousness based on Corine Land Cover of the Mediterranean arc of France:

.. figure:: ../../_static/clc_imperviousness.svg
    :align: center
    :width: 800

On the following picture, we show the discharge results of two uniform calibrations
- green dashed line taking account the imperviousness coefficients, blue dashed line without - 
against observed discharge - orange line.

.. code-block:: python

    code = model.mesh.code[0]
    plt.plot(model.response_data.q[0, :], color='orange', label="Observed discharge")
    plt.plot(model.response.q[0, :], ls='--', color='b', label="Simulated discharge");
    plt.plot(model_imperviousness.response.q[0, :], ls='--', color="g", \
        label="Simulated discharge with imperviousness")
    plt.grid(ls="--", alpha=.7);
    plt.xlabel("Time step (s)");
    plt.ylabel("Discharge ($m^3/s$)")
    plt.title(
        f"Discharges at gauge {code}"
    )
    plt.legend()
    plt.show()