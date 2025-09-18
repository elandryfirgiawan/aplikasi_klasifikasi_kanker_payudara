const background = document.querySelector('.agtc-background');
const letters = ['A', 'G', 'T', 'C'];

function createLetter() {
    const span = document.createElement('span');
    span.classList.add('agtc-letter');
    span.innerText = letters[Math.floor(Math.random() * letters.length)];
    span.style.left = Math.random() * 100 + 'vw';
    span.style.animationDuration = (2 + Math.random() * 3) + 's';
    span.style.fontSize = (12 + Math.random() * 18) + 'px';

    background.appendChild(span);

    setTimeout(() => {
        span.remove();
    }, 5000);
}

// Buat huruf setiap 100ms
setInterval(createLetter, 100);
