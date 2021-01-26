### Tools for the Zwicky Transient Facility (ZTF) Census of the Local Universe (CLU) Experiment
___


![Aitoff Projection of Discovered ZTF CLU Transients](https://github.com/AndyTza/andytza.github.io/blob/master/images/CLU_Map.png?raw=true)


Before running the ZTF CLU scripts, you will need to add your fritz.science credentials. To start using the
fritz.scince API you will need to add to `/scripts/user_info.json` your last name and Command-Line Authentication Token:

```
{
  "user":
    {"user_last_name": "ADD_USER_LAST_NAME_HERE",
     "FritzToken": "ADD_FRITZ_TOKEN_HERE"
    }
}
```

#### Running the `ZTF CLU Luminous Annotation`
____
In ZTF Phase II the ZTF CLU experiment will be a volume-luminosity-limited supernovae sample in hope to find some of the most subluminous transients in the local universe.

This annotator will fetch the photometry of all sources saved to the Caltech ZTF CLU program (`program_id:43`) after the date indicated by the `-d` flag in format `YYYY-MM-DD`. The annotator will go through each source and annotate if they're luminous based on a specific luminosity cut (docs of `luminous_CLU_annotation.py` for more details). To run the ZTF CLU luminosity annotator, see the example below:

```
./luminous_CLU_annotation.py -d 2021-01-01
```
