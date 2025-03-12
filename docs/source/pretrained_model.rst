.. _pretrained_model:

Pretrained Model
=====

- Default pretrained RNA2seg model is available  with the ```pretrained_model``` args in the ```RNA2seg``` class.
    The default pretrained model is trained on the CosMx and merfish datasets accross seven different human organs : breast, colon, lung, melanoma, uterine, prostate and ovarian.
    .. code-block:: console


    rna2seg = RNA2seg(
        device,
        net='unet',
        flow_threshold = 0.9,
        cellbound_flow_threshold = 0.4,
        pretrained_model = "default_pretrained"
    )

- A fine-tune  RNA2seg model on MERFISH hamster brain data can be load with  the ```brain_hamster``` args in the ```RNA2seg``` class.

New train model can be added :) Please open an issue on the github repository.   https://github.com/fish-quant/rna2seg/issues

For other organs or different experimental protocols,  RNA2seg can be fine-tuned on your own data.
See :ref:`userguide` for more details.


