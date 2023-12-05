// Accessibility changes (these are a bit hacky / repetitive while I debug javascript)


// Persistent variable - defaults to pretty font

if (localStorage.getItem('legibility') != 'legible') {
    localStorage.setItem('legibility', 'pretty');
}

// Persistent variable - scale defaults to 1.0

if (localStorage.getItem('fontScale') == null) {
    localStorage.setItem('fontScale', 1.0);
}

// Action this on page load, will use persistent values

// Font scale factor

fontScalingFactor = localStorage.getItem('fontScale')
document.documentElement.style.setProperty('--pst-font-scale-factor', fontScalingFactor);

// Font families

savedStyle = localStorage.getItem('legibility')
if (savedStyle == 'legible') 
    legibleFontSetter('legible')
else 
    legibleFontSetter('pretty')


// Change font size

function fontScaler(scale){
    if(scale == 0.0){
        if (localStorage.getItem('legibility') == 'legible')
            fontScalingFactor = 1.25;
        else
            fontScalingFactor = 1.0;
    }
    else {
        fontScalingFactor = localStorage.getItem('fontScale')
        fontScalingFactor *= scale
    }
    localStorage.setItem('fontScale', fontScalingFactor);
    document.documentElement.style.setProperty('--pst-font-scale-factor', fontScalingFactor);
}

// Change the font family 

function legibleFontSetter(fontType){
    if (fontType == 'legible'){
        legibleFontFamily = getComputedStyle(document.documentElement).getPropertyValue('--pst-font-family-legible');
        legibleFontFamilyH = getComputedStyle(document.documentElement).getPropertyValue('--pst-font-family-legible-headers');
        document.documentElement.style.setProperty('--pst-font-family-base', legibleFontFamily);
        document.documentElement.style.setProperty('--pst-font-family-heading', legibleFontFamilyH);
    }
    else {
        prettyFontFamily = getComputedStyle(document.documentElement).getPropertyValue('--pst-font-family-pretty');
        prettyFontFamilyH = getComputedStyle(document.documentElement).getPropertyValue('--pst-font-family-pretty-headers');
        document.documentElement.style.setProperty('--pst-font-family-base', prettyFontFamily);
        document.documentElement.style.setProperty('--pst-font-family-heading', prettyFontFamilyH);
    }

}

// Toggle between legible and pretty fonts

function legibleFontSwitcher() {
    savedStyle = localStorage.getItem('legibility')
    if (savedStyle == 'pretty') {
        thisStyle = 'legible'
        fontScaler(1.25);
        legibleFontSetter('legible')
    }
    else {
        thisStyle = 'pretty'; 
        fontScaler(0.8);
        legibleFontSetter('pretty')
    }
    localStorage.setItem('legibility', thisStyle);
}


